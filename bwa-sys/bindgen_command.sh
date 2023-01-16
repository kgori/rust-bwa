~/.cargo/bin/bindgen \
  --no-doc-comments \
  --allowlist-function mem_align1_core \
  --allowlist-function mem_sam_pe \
  --allowlist-function mem_opt_init \
  --allowlist-function bwa_idx_load \
  --allowlist-function bwa_idx_destroy \
  --allowlist-function mem_process_seq_pe \
  --allowlist-function mem_process_seqs \
  --allowlist-function mem_align1 \
  --allowlist-function bwa_fill_scmat \
  --allowlist-var "BWA_IDX_.*" \
  --allowlist-var "MEM_F_.*" \
  wrapper.h \
  -o src/lib.rs
