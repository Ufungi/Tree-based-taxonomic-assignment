#!/bin/bash

source ~/anaconda3/etc/profile.d/conda.sh

set -e  # Exit on error

# Set the path to the config.yaml file as a global variable
export CONFIG_GLOBAL="/data/genome/run/snyoo/00_Projects/TM_microbiome_script/Tree-based-taxonomic-assignment/config.yaml"

# Export as environment variable
export WITCH UDANCE FASTROOT

# Load parameters from config file
echo "Loading parameters from config.yaml"
eval $(yq -r 'to_entries | .[] | .key + "=\"" + (.value | tostring) + "\""' ${CONFIG_GLOBAL} | sed -E 's/\(([^)]+)\)/$\1/g')

# Create output directories
mkdir -p ${out}
mkdir -p ${out}/backbone
mkdir -p ${out}/backbone/msa
mkdir -p ${out}/backbone/tree
mkdir -p ${out}/place
mkdir -p ${out}/refine
mkdir -p ${out}/cluster
mkdir -p ${out}/taxonomy

# Set log file
echo "Pipeline started at $(date)" > ${log_file}

# Validate input data (IDs must match exactly)
fasta_ids=$(grep '^>' "${db_fasta}" | sed 's/^>//;s/[[:space:]]//g' | sort | uniq)
tax_ids=$(cut -f1 "${in}/taxonomy.txt" | sed 's/[[:space:]]//g' | sort | uniq)

missing_in_tax=$(comm -23 <(echo "$fasta_ids") <(echo "$tax_ids"))
missing_in_fasta=$(comm -13 <(echo "$fasta_ids") <(echo "$tax_ids"))

if [ -n "$missing_in_tax" ] || [ -n "$missing_in_fasta" ]; then
    [ -n "$missing_in_tax" ] && echo "IDs present in fasta but missing in taxonomy.txt:" && echo "$missing_in_tax"
    [ -n "$missing_in_fasta" ] && echo "IDs present in taxonomy.txt but missing in fasta:" && echo "$missing_in_fasta"
    exit 1
else
    echo "fasta headers and taxonomy.txt first column match exactly."
fi

align_and_build_tree() {
    if [ ! -f ${db_msa} ]; then

        # Align the database sequences
        echo "Running WITCH alignment..." | tee -a ${log_file}
        python3 ${witch} \
        -i ${db_fasta} \
        -d ${out}/backbone/msa \
        -o ${db_prefix}_aligned.fasta \
        -t ${threads}
    else
        echo "Database MSA file already exists: ${db_msa}" | tee -a ${log_file}
    fi

    if [ ! -f ${best_tree} ]; then
        echo "Running VeryFastTree tree inference..." | tee -a ${log_file}

        # Infer a draft tree
        if [ ! -f ${draft_tree} ]; then
            VeryFastTree \
            -nt -gamma -gtr \
            -threads ${threads} \
            ${db_msa} > ${draft_tree}
        else
            echo "Draft tree already exists: ${draft_tree}" | tee -a ${log_file}
        fi

        # Resolve polytomies if resolved_tree does not exist
        if [ ! -f ${resolved_tree} ]; then
            python ${scripts}/resolve_polytomies.py ${draft_tree} ${resolved_tree}
        else
            echo "Resolved tree already exists: ${resolved_tree}" | tee -a ${log_file}
        fi

        # Recalculate branch length and reroot the tree
        if [ ! -f ${best_tree} ]; then
            echo "Running RAxML-NG evaluation..." | tee -a ${log_file}
            raxml-ng \
            --prefix ${wdr}/${organism}/out/backbone/tree/${db_prefix}_aligned.fasta \
            --evaluate --brlen scaled --model GTR --redo \
            --msa ${db_msa} \
            --tree ${resolved_tree} \
            -threads ${threads} $(if [ -n "${outgroup}" ]; then echo "--outgroup ${outgroup}"; fi)
        else
            echo "Best tree already exists: ${best_tree}" | tee -a ${log_file}
        fi
    else
        echo "Best tree file already exists: ${best_tree}" | tee -a ${log_file}
    fi
}

phylogenetic_placement() {

    # Check for special characters (including underscores) in query_fasta
    if grep -q '[^a-zA-Z0-9._/-]' <<< "$(basename ${query_fasta})"; then
        echo "Error: query_fasta filename contains special characters (excluding ., _, -, and /)." | tee -a ${log_file}
        exit 1
    fi

# Chunkify the query sequences

# Calculate chunk_size if set to 'auto' 
if [[ "$chunk_size" == "auto" ]]; then
    n=$(grep -c "^>" "${db_fasta}")
    if (( n / 5 > 500 )); then
        chunk_size=500
    else
        chunk_size=$(( n / 5 ))
    fi
    export chunk_size
    echo "chunk_size set to $chunk_size"
fi
    chunk_exists=$(ls ${out}/place/chunk_*.fasta 2>/dev/null | wc -l)
    if [ "$chunk_exists" -gt 0 ]; then
        echo "Chunk files already exist. Skipping chunkify step." | tee -a ${log_file}
    else
        echo "Chunkifying and placing query sequences..." | tee -a ${log_file}
        gappa prepare chunkify \
        --fasta-path ${query_fasta} \
        --chunk-size ${chunk_size} \
        --chunks-out-dir ${out}/place \
        --abundances-out-dir ${out}/place \
        --threads ${threads}
    fi

    # Place the query sequences
    for fasta_file in ${out}/place/*.fasta; do
        out_jplace=${fasta_file%.fasta}.jplace
        if [ ! -f ${out_jplace} ]; then
            echo "Processing ${fasta_file} with App-SpaM..." | tee -a ${log_file}
            appspam \
            -m SPAMX \
            -p 5 \
            -o ${out_jplace} \
            -t ${best_tree} \
            -s ${db_msa} \
            -q ${fasta_file} \
            --threads ${threads}
        else
            echo "Jplace file already exists: ${out_jplace}" | tee -a ${log_file}
        fi
    done

    # Check query number and add round brackets to the edge of the newick tree
    # See: https://github.com/matthiasblanke/App-SpaM/issues/4
    for out_jplace in ${out}/place/*.jplace; do
        count=$(grep '"p"' ${out_jplace} | wc -w)
        if [ ${count} -eq ${chunk_size} ]; then
            echo "${out_jplace}: Count is ${count}"
        else
            echo "${out_jplace}: Count is not ${chunk_size}"
        fi
        sed 's/"tree":"\(.*\);"/"tree":"(\1);"/g' ${out_jplace} > ${out_jplace%.jplace}.bracket.jplace
    done

    # Unchunkify the query sequences (output = query.jplace)
    echo "Unchunkifying query sequences..." | tee -a ${log_file}
    gappa prepare unchunkify \
    --abundances-path ${out}/place/abundances_query.json \
    --jplace-path ${out}/place/*.bracket.jplace \
    --allow-file-overwriting \
    --threads ${threads} \
    --out-dir ${out}/place
}

fix_jplace_file() {
    echo "Examining query.jplace and fixing any issues..." | tee -a ${log_file}

    # Try to examine original file first and fix if needed
    if ! gappa examine info --jplace-path ${out}/place/query.jplace; then
        echo "Removing placements with nan values using sed..." | tee -a ${log_file}

        # Remove only the specific "nan" entries in the "p" field
        sed -e '/"p":\s*\[\s*\[[^]]*nan[^]]*\]\s*\]/d' \
            ${out}/place/query.jplace > ${out}/place/query.clean.jplace

        # Verify the fixed file
        if ! gappa examine info --jplace-path ${out}/place/query.clean.jplace; then
            echo "Error persists in query.clean.jplace. Manual intervention required." | tee -a ${log_file}
            exit 1
        fi

        echo "File query.clean.jplace has been fixed and is ready for processing." | tee -a ${log_file}
    else
        echo "File query.jplace passed validation and is ready for processing." | tee -a ${log_file}
    fi

    # Define final jplace file to use
    if [ -f "${out}/place/query.clean.jplace" ]; then
        export fixed_jplace="${out}/place/query.clean.jplace"
    else
        export fixed_jplace="${out}/place/query.jplace"
    fi

echo "Split multiplicity..." | tee -a ${log_file}
    
# split multiplicity
python3 - <<EOF
import json

input_file = "${fixed_jplace}"
output_file = "${out}/place/query.split.jplace"

with open(input_file) as f:
    data = json.load(f)

new_placements = []
for pl in data["placements"]:
    if len(pl["n"]) > 1:
        for n in pl["n"]:
            new_placements.append({"p": pl["p"], "n": [n]})
    else:
        new_placements.append(pl)

data["placements"] = new_placements

with open(output_file, "w") as f:
    json.dump(data, f, indent=2)
EOF

}

tree_refinement() {
    echo "Preparing and running UDance analysis..." | tee -a ${log_file}

    # Create UDance directories
    mkdir -p ${out}/refine/{alignments,output/{placement,trimmed}}

    # Prepare sequences
    echo "Preparing sequences for UDance..." | tee -a ${log_file}

    # Pad sequences and prepare files
    # uDance originally requires aligned sequences, but can accept unaligned sequences with equal length
    python ${scripts}/seq_length_sync.py -i ${query_fasta} -o ${out}/refine/output/placement/query.fa
    
    # Combine db and query fasta and verify the line count
    cat ${db_fasta} ${query_fasta} > ${out}/refine/combined.fasta
    if [ $(wc -l < ${out}/refine/combined.fasta) -ne $(( $(wc -l < ${db_fasta}) + $(wc -l < ${query_fasta}) )) ]; then
        echo "Error: File concatenation failed - line count mismatch" | tee -a ${log_file}
        exit 1
    fi

    python ${scripts}/seq_length_sync.py -i ${out}/refine/combined.fasta -o ${out}/refine/alignments/backbone.fa
    rm ${out}/refine/combined.fasta
    cp ${out}/refine/alignments/backbone.fa ${out}/refine/output/placement/backbone.fa
    cp ${out}/refine/alignments/backbone.fa ${out}/refine/output/trimmed/backbone.fa

    cp ${best_tree} ${out}/refine/output/backbone.nwk

    cp ${fixed_jplace} ${out}/refine/output/placement.jplace
    # Replace all occurrences of "nm" with "n" in placement.jplace
    # refer to decompose.py line 25
    sed -i 's/"nm":/"n":/g' ${out}/refine/output/placement.jplace

    # Activate uDance environment
    echo "Running UDance analysis..." | tee -a ${log_file}
    conda activate udance

    # Change to uDance directory
    pushd $udance > /dev/null

    # Update workdir in udance config file
    sed -i 's|workdir:.*|workdir: "'${out}/refine/'"|' config.yaml

    # Run uDance
    snakemake -c ${threads} --configfile config.yaml --snakefile udance.smk all --unlock
    snakemake -c ${threads} --configfile config.yaml --snakefile udance.smk all

    # Return to the original directory
    popd > /dev/null
}

cluster_analysis() {
    echo "Starting cluster analysis..." | tee -a ${log_file}

    # Get the list of subfolders matching the pattern
    subfolders=$(find ${out}/refine/output/udance -mindepth 1 -maxdepth 1 -type d -regex '.*/[0-9]+')

    # Process each subfolder
    for folder in ${subfolders}; do
        echo "Processing folder: ${folder}" | tee -a ${log_file}
        cd "${folder}"

        # Check if the required tree file exists
        if [ ! -f ./backbone/raxml.expanded.nwk ]; then
            echo "File ./backbone/raxml.expanded.nwk not found in ${folder}. Skipping." | tee -a ${log_file}
            cd - > /dev/null
            continue
        fi

        # Root the tree using FastRoot
        python $fastroot \
        -i ./backbone/raxml.expanded.nwk \
        -m MP \
        -o ./raxml.expanded.nwk.rooted

        # Perform clustering using TreeCluster
        cluster_file="${out}/cluster/cluster_$(basename ${folder}).tsv"
        TreeCluster.py -i ./raxml.expanded.nwk.rooted -t 0.02 -m single_linkage > ${cluster_file}

        # Copy the cluster file to the output directory
        aln_file="${out}/cluster/shrunk.aln_$(basename ${folder}).fasta"
        cp ./backbone/shrunk.aln ${aln_file}

        # Copy the tree file to the output directory
        tree_file="${out}/cluster/raxml.expanded.nwk.rooted_$(basename ${folder}).nwk"
        cp ./raxml.expanded.nwk.rooted ${tree_file}

        cd - > /dev/null
    done
}

detect_long_branches() {
    echo "Detecting long-branched tips..." | tee -a ${log_file}

    # Change to the cluster output directory
    cd ${out}/cluster/

    # Process alignment files
    for aln_file in shrunk.aln_*.fasta; do
        suffix=${aln_file#shrunk.aln_}
        suffix=${suffix%.fasta}
        tree_file="raxml.expanded.nwk.rooted_${suffix}.nwk"

        if [[ -f ${tree_file} ]]; then
            echo "Running TreeShrink on ${aln_file}..." | tee -a ${log_file}
            run_treeshrink.py -a ${aln_file} -t ${tree_file} -o ./ -O ./shrunk_${suffix} -q 0.25 -b 5 -f -c
        else
            echo "Tree file ${tree_file} not found for alignment file ${aln_file}. Skipping." | tee -a ${log_file}
        fi
    done

    # Define the output TSV file
    long="long-branched_tips.tsv"
    > ${long}

    # Process TreeShrink output files
    for txt_file in shrunk_*.txt; do
        if [[ ${txt_file} =~ ^shrunk_[^_]+\.txt$ ]]; then
            suffix=${txt_file#shrunk_}
            suffix=${suffix%.txt}
            while IFS=$'\t' read -ra tips; do
                for tip in "${tips[@]}"; do
                    echo -e "${suffix}\t${tip}" >> ${long}
                done
            done < ${txt_file}
        else
            echo "Skipping file ${txt_file} as it does not match the strict pattern." | tee -a ${log_file}
        fi
    done
}

run_taxonomy_script() {
    echo "Running LCA taxonomy script..." | tee -a ${log_file}
    Rscript ${scripts}/LCA_taxonomy.R --in ${in} --path ${out} --DB ${db_fasta}
}

print_help() {
    cat <<EOF
[Pipeline Usage]

Usage: $0 [step|range|all|-h|--help]

Options:
  step      : Run a specific step only (e.g., 3)
  range     : Run a range of steps (e.g., 3-6)
  all       : Run all steps in order (default if no argument)
  -h, --help: Show this help message

[Step List & Description]
  1: align_and_build_tree
     - Aligns the database sequences and builds the initial phylogenetic tree.
  2: phylogenetic_placement
     - Places query sequences onto the reference tree (e.g., using App-SpaM).
  3: fix_jplace_file
     - Checks and fixes errors in the jplace file.
  4: tree_refinement
     - Refines the tree and performs additional analysis with UDance.
  5: cluster_analysis
     - Performs tree-based clustering and generates cluster files.
  6: detect_long_branches
     - Detects long-branched tips using TreeShrink and records them.
  7: run_taxonomy_script
     - Runs the R script for LCA-based taxonomy assignment and final output.

[Examples]
  $0             # Run all steps in order
  $0 all         # Run all steps in order
  $0 3           # Run only step 3
  $0 3-6         # Run steps 3, 4, 5, and 6 in order
  $0 -h          # Show this help message

[Notes]
- Step numbers range from 1 to 7.
- If you specify a range, the script will run all steps from the start to the end number (inclusive).
- If a step fails, the pipeline continues to the next step.
- Parameters are automatically loaded from config.yaml.
- Logs are written to \$log_file.

EOF
}

main() {
    set -e  # Exit immediately if a command exits with a non-zero status

    steps=(align_and_build_tree phylogenetic_placement fix_jplace_file tree_refinement cluster_analysis detect_long_branches run_taxonomy_script)

    # Help
    if [[ "$1" == "-h" || "$1" == "--help" ]]; then
        print_help
        exit 0
    fi

    # Run all steps (all or no argument)
    if [[ -z "$1" || "$1" == "all" ]]; then
        for i in "${!steps[@]}"; do
            ${steps[$i]}
        done
        echo "Pipeline completed successfully at $(date)" | tee -a ${log_file}
        return
    fi

    # Range (e.g., 3-6)
    if [[ "$1" =~ ^([0-9]+)-([0-9]+)$ ]]; then
        start=${BASH_REMATCH[1]}
        end=${BASH_REMATCH[2]}
        for ((i=start; i<=end; i++)); do
            idx=$((i-1))
            ${steps[$idx]}
        done
        echo "Pipeline completed successfully at $(date)" | tee -a ${log_file}
        return
    fi

    # Single step (e.g., 3)
    if [[ "$1" =~ ^[0-9]+$ ]]; then
        idx=$(( $1 - 1 ))
        ${steps[$idx]}
        echo "Pipeline completed successfully at $(date)" | tee -a ${log_file}
        return
    fi

    # Otherwise, print help
    print_help
    exit 1
}

main "$1"
