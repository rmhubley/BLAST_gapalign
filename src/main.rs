//! # BLAST_GapAlign
//!
//! This is a Rust port of the BLAST ALIGN_EX function for gapped alignment 
//! extension. It is based on the original C code found in the 2.14.1 release 
//! of NCBI BLAST with minimal modifications to simplify the data structures.
//! This is **not** a full port of the BLAST algorithm, which involves many more
//! components (seed finding, ungapped extension, statistics generation, sequence
//! management etc). However, it provides the same means of generating a greedy
//! alignment of two subsequences given an anchor (seed), a scoring system 
//! (matrix, and affine gap penalties) and a x-drop threshold for termining the 
//! extension of the alignment.
//!
//! ## Features
//! - A simple API for greedy seeded alignment extension.
//! - Support for custom scoring matrices and affine gap penalties.
//!
//! ## Example
//!
//!
//! ## References
//!
//!  Original C code: <https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.14.1/ncbi-blast-2.14.1+-src.tar.gz>
//!  file: ncbi-blast-2.14.1+-src/c++/src/algo/blast/core/blast_gapalign.c
//!
//!  Zhang, Zheng, et al. 
//!  "*A greedy algorithm for aligning DNA sequences.*" 
//!  Journal of Computational biology 7.1-2 (2000): 203-214.
//!
//!  ## TODO
//!  - Modify the align_ex function to support extensions on the reverse strand of
//!    sequence_b.
//!  - Document the align_ex algorithm in more detail.
//!
use std::cmp::{max, min};
use std::i32::MIN;

/// Converts a DNA character to contigous numeric encoding.
fn dna_base_to_index(base: u8) -> usize {
    match base {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        b'N' => 4,
        _ => panic!("Invalid DNA base: {}", base),
    }
}

// Edit operations ( for traceback )
pub const OP_SUBSTITUTION:u8  = 0x03;
pub const OP_GAP_IN_A:u8       = 0x00; 
pub const OP_GAP_IN_B:u8    = 0x06;
pub const OP_EXTEND_GAP_A:u8 = 0x10;
pub const OP_EXTEND_GAP_B:u8  = 0x40;

/// Scoring parameters for gapped alignment
/// - A slightly modified version of the original NCBI BLAST scoring parameters
///   structure. For this simple API it makes more sense to combine the scoring
///   matrix, gap penalties, and x-drop threshold into a single structure.
struct BlastScoringParameters {
    scoring_matrix: Vec<Vec<i32>>,       // Standard scoring matrix
    gap_open: i32,
    gap_extend: i32,
    gap_x_dropoff: i32,
}

// Individual cell of the dynamic programming matrix
#[derive(Clone)]
struct BlastGapDP {
    best: i32,        // Best score for the cell
    best_gap: i32,    // Best score for a gap (insertion/deletion)
}

// Dynamic programming matrix for gapped alignment
struct BlastGapAlignStruct {
    dp_mem: Vec<BlastGapDP>, // Holds a vector of BlastGapDP cells
    dp_mem_alloc: usize,     // Size of the allocated memory for dp_mem
}

// Print a row/diagonal/anti-diagonal of the DP matrix -- not sure yet
#[allow(dead_code)]
fn print_scores_range(scores: &Vec<BlastGapDP>, start: usize, end: usize) {
    // Ensure the range is valid
    if start >= scores.len() || end > scores.len() || start >= end {
        println!("Invalid range: start = {}, end = {}, length = {}", start, end, scores.len());
        return;
    }

    // Print the "best" field
    print!("best: ");
    for i in start..end {
        print!("{}", scores[i].best);
        if i < end - 1 {
            print!(", ");
        }
    }
    println!();

    // Print the "best_gap" field
    print!("best_gap: ");
    for i in start..end {
        print!("{}", scores[i].best_gap);
        if i < end - 1 {
            print!(", ");
        }
    }
    println!();
}

//
// A Rust port of the NCBI BLAST blast_gapalign.c:ALIGN_EX function
//
// Notes on conventions:
//
// DP Matrix orientation:
//        seq A
//       +---------
//     s |
//     e |
//     q |
//       |
//     B |
//       
fn align_ex(
    a: &[u8],                                // Sequence A
    b: &[u8],                                // Sequence B
    m: i32,                                  // Length of sequence A 
    n: i32,                                  // Length of sequence B
    a_offset: &mut i32,                      // Best alignment endpoint in A (output)
    b_offset: &mut i32,                      // Best alignment endpoint in B (output)
    gap_align: &mut BlastGapAlignStruct,     // Structure containing alignment state
    score_params: &BlastScoringParameters,   // Parameters for scoring (gap penalties, etc.)
    reverse_sequence: bool,                  // Reverse the sequence direction (e.g. extend left)
) -> (i32, Vec<u8>) {                        // Returns best alignment score, and the edit script

    if n <= 0 || m <= 0 {
        return (0, Vec::new());
    }

    if n >= b.len() as i32 || m >= a.len() as i32 {
        panic!("align_ex: Sequence length exceeds buffer size");
    }

//println!("Aligning sequences of length {} and {}", m, n);
//println!("reverse_sequence: {}", reverse_sequence);
//println!("gap_open: {}, gap_extend: {}", score_params.gap_open, score_params.gap_extend);

    if score_params.gap_open <= 0 {
        panic!("align_ex: Gap open penalty must be positive integer greater than 0.");
    }

    if score_params.gap_extend <= 0 {
        panic!("align_ex: Gap extension penalty must be positive integer greater than 0.");
    }

    // --- Step 1: Initialization ---
    *a_offset = 0;
    *b_offset = 0;

    let mut x_dropoff = score_params.gap_x_dropoff;
    let gap_open = score_params.gap_open;
    let gap_extend = score_params.gap_extend;

    let gap_open_extend = gap_open + gap_extend;

    if x_dropoff < gap_open_extend {
        x_dropoff = gap_open_extend;
    }

    // Allocate initial memory for traceback rows and start offsets
    let mut edit_script: Vec<Vec<u8>> = Vec::new();
    let mut edit_start_offset: Vec<usize> = Vec::new();

   //
    // The smaller the extension penalty the longer a dp row will need to be
    // To be conservative we will allocate enough memory for the worst case
    //   - all gaps in in B until the xdrop is reached = xdrop/gap_extend
    //     + 3 for safety?
    let num_extra_cells = ((x_dropoff / score_params.gap_extend) + 3) as usize;
//println!("num_extra_cells: {}", num_extra_cells);

    // Allocate memory for DP matrix
    //  It appears that the minimum allocation size is 100
    if num_extra_cells > gap_align.dp_mem_alloc {
        // Not sure why the +100 is necessary
        gap_align.dp_mem_alloc = num_extra_cells + 100; 
        gap_align.dp_mem.resize(
            gap_align.dp_mem_alloc,
            BlastGapDP {
                best: MIN,
                best_gap: MIN,
            },
        );
    }

    // Setup pointers to DP row and edit_script row
    edit_script.push(vec![0;num_extra_cells]);
    edit_start_offset.push(0);
    let score_array = &mut gap_align.dp_mem;

    // Initialize the first cell of the dp_row
    let mut dp_row_used = 1;
    let mut score = -gap_open_extend;
    score_array[0].best = 0;
    score_array[0].best_gap = -gap_open_extend;

    // Initialize the first row of the DP matrix
    // NOTE: this used to be "1..=n" translated from the original C code
    //       but this seems dangerous as we do know the real length of the
    //       score_array and n could be larger than it.  Also, dp_row_used 
    //       should also not be larger than then n.  So we should use the 
    //       min of the two, which means we may not reach the x_dropoff in this loop.
    for i in 1..min(n,gap_align.dp_mem_alloc as i32) as usize {
//println!("  setting sc_used: i: {}, score: {}", i, score);
        if score < -x_dropoff {
            break;
        }
        score_array[i].best = score;
        score_array[i].best_gap = score - gap_open_extend;
        score -= gap_extend;
        dp_row_used += 1;
    }
//println!("First row of DP matrix[{}..] (allocated={}, used={}): ", 0, gap_align.dp_mem_alloc, dp_row_used);
//print_scores_range(&score_array, 0, dp_row_used);

    // So b_size is the number of cells in the dp_row that are in-use
    let mut b_size = dp_row_used as usize;
    let mut best_score = 0;
    let mut first_b_index = 0;

    // --- Step 2: Outer loop over sequence A ---
    for a_index in 1..=m as usize {
//println!("Starting a_index: {}", a_index);
        // TODO: This doesn't need to be this large...
        edit_script.push(vec![0;b_size+3]);
        edit_start_offset.push(first_b_index);

        let edit_script_row = &mut edit_script[a_index];
        let mut last_b_index = first_b_index;

        // Reset score to minimum int value
        score = MIN;
        let mut score_gap_row = MIN;

//println!("  b_size: {}, first_b_index: {}, last_b_index: {}, score: {}, score_gap_row: {}", b_size, first_b_index, last_b_index, score, score_gap_row);

        // --- Step 3: Inner loop over sequence B ---
        for b_index in first_b_index..b_size {
//println!("  Starting b_index: {} and score: {}",b_index,score);

            let matrix_score = if reverse_sequence { 
//println!("--> n={}, b_index={}", n, b_index);
                if  n - (b_index as i32) == 0  {
                    MIN
                } else {
                    score_params.scoring_matrix[dna_base_to_index(a[(m - (a_index as i32)) as usize])]
                                        [dna_base_to_index(b[(n - (b_index as i32) - 1) as usize])]
                }
           } else {
                if b_index == n as usize {
                    MIN
                } else {
                    score_params.scoring_matrix[dna_base_to_index(a[a_index])]
                                        [dna_base_to_index(b[b_index + 1])]
                }
            };

            let next_score = score_array[b_index].best.saturating_add(matrix_score);

            // Previous Insertion score 
            let mut score_gap_col = score_array[b_index].best_gap;

/*
if reverse_sequence {
    let a_idx = (m - (a_index as i32)) as usize;
    let a_base = a[a_idx];
    let b_idx = if n - (b_index as i32) == 0 {
        MIN
    } else {
        (n - (b_index as i32) - 1)
    };
    let b_base = if n - (b_index as i32) == 0 {
        25 // '%'
    }else {
         b[b_idx as usize]
    };
    println!("    reverse_sequence: m-a_index: {} = {}, n-b_index: {} = {} matrix_score = {}, next_score = {}",
             a_idx, a_base as char, b_idx, b_base as char, matrix_score, next_score);
}else {
    let a_base = a[a_index];
    let b_base = if b_index == n as usize {
        25 // '%' 
    } else {
        b[b_index+1]
    };
    println!("    reverse_sequence: a_index: {} = {}, b_index+1: {} = {} matrix_score = {}, next_score = {}",
             a_index, a_base as char, b_index+1, b_base as char, matrix_score, next_score);
}
*/

            let mut script = OP_SUBSTITUTION;
//println!("    ******* RIGHT BEFORE score = {} score_gap_col = {} and score_gap_row = {}\n", score, score_gap_col, score_gap_row);
            if score < score_gap_col {
//println!("    ====score {} < score_gap_col {} setting to that", score, score_gap_col);
                score = score_gap_col;
                script = OP_GAP_IN_B;
            }
            if score < score_gap_row {
//println!("    ====score {} < score_gap_row {} setting to that", score, score_gap_col);
                score = score_gap_row;
                script = OP_GAP_IN_A;
            }
//println!("1: best_score = {}, score = {}, score_gap_col = {}", best_score, score, score_gap_col);

//println!("    bestscore: {}, score: {} (diff = {}), score_gap_col = {}, score_gap_row = {}", best_score, score, best_score - score, score_gap_col, score_gap_row);
            if best_score - score > x_dropoff {
                if b_index == first_b_index {
                    first_b_index += 1;
                } else {
                    score_array[b_index].best = MIN;
                }
            } else {
                last_b_index = b_index;
//println!("    setting last_b_index: {}", last_b_index);
                if score > best_score {
                    best_score = score;
                    *a_offset = a_index as i32;
                    *b_offset = b_index as i32;
                }
     
//println!("B:  score_gap_col = {}, score_gap_row = {}, gap_extend = {}", score_gap_col, score_gap_row, gap_extend);
                score_gap_col = score_gap_col.saturating_sub(gap_extend);
                if score_gap_col < score - gap_open_extend {
                    //println!("Set best_gap to {} (gap_open_extend = {}", score - gap_open_extend, gap_open_extend);
                    score_array[b_index].best_gap = score - gap_open_extend;
                } else {
                    //println!("Set best_gap to {}", score_gap_col);
                    score_array[b_index].best_gap = score_gap_col;
                    script |= OP_EXTEND_GAP_B;
                }
                score_gap_row = score_gap_row.saturating_sub(gap_extend);
//println!("    ####### score_gap_row {} <? score {} - gap_open_extend {} = {}", score_gap_row, score, gap_open_extend, score - gap_open_extend);
                if score_gap_row < score - gap_open_extend {
//println!("    ####### YES");
                    score_gap_row = score - gap_open_extend;  // Open a new deletion
                } else {
                    script |= OP_EXTEND_GAP_A;
                }
                score_array[b_index].best = score;
            }
//println!("A:  score_gap_col = {}, score_gap_row = {}, gap_extend = {}", score_gap_col, score_gap_row, gap_extend);
//println!("    Setting score = {}", next_score);
            score = next_score;

            // Save the traceback script
//println!("    Setting edit script row b_index={} first_b_index={} script={}", b_index, first_b_index, script);
            edit_script_row[b_index] = script as u8;

        } // b-loop

        // TODO Revist this logic along with fence_hit
        if first_b_index == b_size {
//println!("breaking because first_b_index == b_size");
            break;
        }

        if last_b_index + num_extra_cells + 3 >= score_array.len() {
            // Dynamically expand score_array
            let new_size = max(b_size + 100, score_array.len() * 2);
//println!("Expanding score_array to size {}", new_size);
            score_array.resize(
                        new_size,
                        BlastGapDP {
                            best: MIN,
                            best_gap: MIN,
                        },
                    );
        }
        // NEW
        //if ( b_size >= edit_script_row.len() ) {
        //            // Dynamically expand edit_script_row
        //            let new_size = max(b_size + 100, edit_script_row.len() * 2);
//println!("Expanding edit_script_row to size {}", new_size);
 //                   edit_script_row.resize(new_size, 0);
 //               }



        // If the row extends beyond the current `b_size`, expand it with gap entries
        if last_b_index < b_size - 1 {
            b_size = last_b_index + 1;
        }else {
//println!("HERE: score_gap_row = {}, best_score = {}, x_dropoff = {}, b_size = {}", score_gap_row, best_score, x_dropoff, b_size);
            while score_gap_row >= (best_score - x_dropoff) && b_size <= n as usize {
                score_array[b_size].best = score_gap_row;
                score_array[b_size].best_gap = score_gap_row - gap_open_extend;
                score_gap_row -= gap_extend;

                if b_size - first_b_index >= edit_script_row.len() {
                    panic!(
                        "Index out of bounds in edit_script_row: b_size = {}, first_b_index = {}, row size = {}",
                        b_size, first_b_index, edit_script_row.len()
                    );
                }
                edit_script_row[b_size] = OP_GAP_IN_A;
                b_size += 1;
//println!("1expanding bsize to {}", b_size);
            }
        }
        if b_size <= n as usize {
            score_array[b_size].best = MIN;
            score_array[b_size].best_gap = MIN;
            b_size += 1;
//println!("2expanding bsize to {}", b_size);
        }

//println!("  Edit Script Row: {:?}", edit_script_row);
//println!("  Edit Script for row {}:", a_index);

//println!("Bsize is now {}", b_size);

/*
// Iterate through the row and print each operation with its corresponding b_index
    for (b, &op) in edit_script_row[first_b_index..b_size].iter().enumerate() {
        let b_index = b + first_b_index;
         match op & 0x07 {
            OP_SUBSTITUTION => println!("    b_index = {}: Substitution/Match == {}", b_index, op),
            OP_GAP_IN_A => { if (op & OP_EXTEND_GAP_A) > 0 {
                               println!("    b_index = {}: Extend Gap in A == {}", b_index, op);
                            } else {
                               println!("    b_index = {}: Gap in A == {}", b_index, op);
                            }
                         },
            OP_GAP_IN_B => { if (op & OP_EXTEND_GAP_B) > 0 {
                               println!("    b_index = {}: Extend Gap in B == {}", b_index, op);
                            } else {
                               println!("    b_index = {}: Gap in B == {}", b_index, op);
                            }
                         },
            _ => println!("    b_index = {}: Unknown Operation == {}", b_index, op),
        };
    }
//println!("Row {} of DP matrix[{}..] (allocated={}, used={}-{}): ", a_index, first_b_index, gap_align.dp_mem_alloc, first_b_index, b_size);
//print_scores_range(&score_array, first_b_index, b_size);
*/

    } // a-loop

    // --- Step 4: Traceback ---
    let mut a_index = *a_offset as usize;
    let mut b_index = *b_offset as usize;
    let mut script = OP_SUBSTITUTION;
    let mut edit_sc = Vec::new();
    while a_index > 0 || b_index > 0 {
//println!("starting a_index: {}, b_index: {}, script: {}", a_index, b_index, script);
        //let next_script = edit_script[a_index-1][b_index];
        let next_script = edit_script[a_index][b_index];

        match script {
            OP_GAP_IN_A => { script = next_script & 0x07; 
                           if (next_script & OP_EXTEND_GAP_A) > 0 {
                               script = OP_GAP_IN_A;
                           }
                         },
            OP_GAP_IN_B => { script = next_script & 0x07; 
                           if (next_script & OP_EXTEND_GAP_B) > 0 {
                               script = OP_GAP_IN_B;
                           }
                         },
            _ => {
                script = next_script & 0x07;
            },
        }

//println!("   -- a_index: {}, b_index: {}, script: {}, next_script: {}", a_index, b_index, script, next_script);
        if script == OP_GAP_IN_A {
            b_index -= 1;
        } else if script == OP_GAP_IN_B {
            a_index -= 1;
        } else {
            a_index -= 1;
            b_index -= 1;
        }

        edit_sc.push(script);
    }

    (
        best_score,
        edit_sc
    )
}

fn reconstruct_alignment(
    sequence_a: &[u8],
    sequence_b: &[u8],
    edit_block: Vec<u8>,
    start_a: usize,
    start_b: usize,
) -> (Vec<u8>, Vec<u8>) {
    let mut aligned_a = Vec::new();
    let mut aligned_b = Vec::new();

    let mut a_idx = start_a;
    let mut b_idx = start_b;

    for op in edit_block {
        match op & 0x07 {
            OP_SUBSTITUTION => {
                aligned_a.push(sequence_a[a_idx]);
                aligned_b.push(sequence_b[b_idx]);
                a_idx += 1;
                b_idx += 1;
            }
            OP_GAP_IN_A => {
                aligned_a.push(b'-');
                aligned_b.push(sequence_b[b_idx]);
                b_idx += 1;
            }
            OP_GAP_IN_B => {
                aligned_a.push(sequence_a[a_idx]);
                aligned_b.push(b'-');
                a_idx += 1;
            }
            _ => (),
        }
    }

    (aligned_a, aligned_b)
}


#[allow(dead_code)]
fn align_sequences_with_seed(
    sequence_a: &[u8],
    sequence_b: &[u8],
    seed_a_offset: usize,
    seed_b_offset: usize,
    gap_open: i32,
    gap_extend: i32,
    scoring_matrix: Vec<Vec<i32>>,
    x_dropoff: i32,
) -> (
    String, // Aligned sequence_a
    String, // Aligned sequence_b
    i32,   // Left alignment score
    i32,   // Right alignment score
    (usize, usize), // Sequence A range (start, end)
    (usize, usize), // Sequence B range (start, end)
) {
    // --- Create the alignment structure ---
    let mut gap_align = BlastGapAlignStruct {
        dp_mem: Vec::new(),
        dp_mem_alloc: 0,
    };

    let scoring_params = BlastScoringParameters {
        scoring_matrix: scoring_matrix,
        gap_open,
        gap_extend,
        gap_x_dropoff: x_dropoff,
    };

    // --- Left alignment extension ---
    let mut a_offset = 0;
    let mut b_offset = 0;

    let (left_score, left_edit_script) = align_ex(
        sequence_a,
        sequence_b,
        seed_a_offset as i32 + 1,
        seed_b_offset as i32 + 1,
        &mut a_offset,
        &mut b_offset,
        &mut gap_align,
        &scoring_params,
        true, // Reverse sequence direction (left extension)
    );

    // The left extension includes the base at the seed location
    // so a/b_offset should always be at least 1
    let a_start = seed_a_offset.saturating_sub((a_offset-1) as usize);
    let b_start = seed_b_offset.saturating_sub((b_offset-1) as usize);

    // Reconstruct left aligned sequences
    let (left_a, left_b) =
        reconstruct_alignment(sequence_a, sequence_b, left_edit_script, a_start as usize, b_start as usize);

    // --- Right alignment extension ---
    a_offset = 0;
    b_offset = 0;

    let (right_score, mut right_edit_script) = align_ex(
        &sequence_a[(seed_a_offset as usize)..],
        &sequence_b[(seed_b_offset as usize)..],
        (sequence_a.len() - (seed_a_offset as usize) - 1) as i32,
        (sequence_b.len() - (seed_b_offset as usize) - 1) as i32,
        &mut a_offset,
        &mut b_offset,
        &mut gap_align,
        &scoring_params,
        false, // Forward sequence direction (right extension)
    );

    let a_end = seed_a_offset + (a_offset as usize);
    let b_end = seed_b_offset + (b_offset as usize);

    // Reconstruct left aligned sequences
    right_edit_script.reverse();
    let (right_a, right_b) = reconstruct_alignment(&sequence_a, &sequence_b, right_edit_script, (seed_a_offset+1) as usize, (seed_b_offset+1) as usize);

    let aligned_a = String::from_utf8(left_a.into_iter().chain(right_a).collect()).unwrap();
    let aligned_b = String::from_utf8(left_b.into_iter().chain(right_b).collect()).unwrap();


    (
        aligned_a,
        aligned_b,
        left_score,
        right_score,
        (a_start, a_end),
        (b_start, b_end),
    )
}



fn main() {
    // --- Step 1: Define sequences ---
    // One change seqA: G to seqB: A
    // rmblast says: 175 comparison matrix (symmetric)
    // rmblast says: 166 for 25p41g-Hpad.matrix
    let sequence_a: Vec<u8> = b"TTCGTACGCGATCGAGTACC".to_vec();
    let sequence_b: Vec<u8> = b"TTCGTACGCAATCGAGTACC".to_vec();
    
    // Seed location
    let seed_a_offset = 15;
    let seed_b_offset = 15; 

    // --- Step 2: Define scoring parameters ---
    let gap_open = 20;
    let gap_extend = 5;

    // Create a 5x5 scoring matrix for DNA bases (A, C, G, T, N) 
    // -- 175
    let scoring_matrix = vec![
      //   A   C   G   T   N
      vec![9, -15, -6, -17, -1], // A
      vec![-15, 10, -15, -6, -1], // C
      vec![-6, -15, 10, -15, -1], // G
      vec![-17, -6, -15, 9, -1], // T
      vec![-1, -1, -1, -1, -1], // N 
    ];

    // So scoring matrix orientation is
    // the same as NCBI matrices
    //      subject (b)
    // q
    // u
    // e
    // r
    // y
    // (a)

    /*
    // tt1.fa tt2.fa
    //25p41g.matrix -- 166
    let scoring_matrix = vec![
      // A   C  G   T   N                                
      vec![   8, -13,  -2, -15,  -1],                             
      vec![ -13,  10, -13,  -6,  -1],                                
      vec![  -6, -13,  10, -13,  -1],                                
      vec![ -15,  -2, -13,   8,  -1],                                
      vec![  -1,  -1   -1,  -1,  -1],
    ];
    */


    // --- Step 3: Create the alignment structure ---
    let mut gap_align = BlastGapAlignStruct {
        dp_mem: Vec::new(),
        dp_mem_alloc: 0,
    };

    let scoring_params = BlastScoringParameters {
        scoring_matrix: scoring_matrix,
        gap_open: gap_open,
        gap_extend: gap_extend,
        gap_x_dropoff: 50, // X-dropoff threshold (xdrop_gap_final in NCBI Blast)
    };

    let mut total_score = 0;

    // --- Step 4: Run the left alignment extension ---
    let mut a_offset = 0;
    let mut b_offset = 0;

    let (left_score, left_edit_script) = align_ex(
        &sequence_a,
        &sequence_b,
        seed_a_offset+1 as i32,
        seed_b_offset+1 as i32,
        &mut a_offset,
        &mut b_offset,
        &mut gap_align,
        &scoring_params,
        true, // reverse sequence direction ( left extension ) 
    );

    // --- Step 5: Print results ---
    println!("Best alignment score: {}", left_score);
    println!("Alignment ends at:");
    println!("  Sequence A offset: {}", a_offset);
    println!("  Sequence B offset: {}", b_offset);
    let a_start = seed_a_offset - a_offset + 1;
    let b_start = seed_b_offset - b_offset + 1;
    total_score += left_score;

    // Reconstruct left aligned sequences
    let (left_a, left_b) = reconstruct_alignment(&sequence_a, &sequence_b, left_edit_script, a_start as usize, b_start as usize);

    /*
    println!("Edit operations (traceback):");
    for (op, length) in &edit_block.ops {
        match op & 0x07 {
            OP_SUBSTITUTION => println!("  Substitution: {}", length),
            OP_GAP_IN_A => println!("  Gap in A (insertion): {}", length),
            OP_GAP_IN_B => println!("  Gap in B (deletion): {}", length),
            OP_EXTEND_GAP_A => println!("  Extend gap in A: {}", length),
            OP_EXTEND_GAP_B => println!("  Extend gap in B: {}", length),
            _ => println!("  Unknown operation: {} ({})", op, length),
        }
    }
    */


    // --- Step 6: Run the right alignment extension ---
    a_offset = 0;
    b_offset = 0;

    println!("#\n# RUNNING RIGHT\n#");

    let (right_best_score, mut right_edit_script) = align_ex(
        &sequence_a[(seed_a_offset as usize)..],
        &sequence_b[(seed_b_offset as usize)..],
        (sequence_a.len() - (seed_a_offset as usize) - 1) as i32,
        (sequence_b.len() - (seed_b_offset as usize) - 1) as i32,
        &mut a_offset,
        &mut b_offset,
        &mut gap_align,
        &scoring_params,
        false, // reverse_sequence
    );

    // --- Step 7: Print results ---
    println!("Best alignment score: {}", right_best_score);
    println!("Alignment ends at:");
    println!("  Sequence A offset: {}", a_offset);
    println!("  Sequence B offset: {}", b_offset);
    let a_end = seed_a_offset + a_offset;
    let b_end = seed_b_offset + b_offset;
    total_score += right_best_score;

/*
    println!("Edit operations (traceback):");
    for (op, length) in &edit_block.ops {
        match op & 0x07 {
            OP_SUBSTITUTION => println!("  Substitution: {}", length),
            OP_GAP_IN_A => println!("  Gap in A (insertion): {}", length),
            OP_GAP_IN_B => println!("  Gap in B (deletion): {}", length),
            OP_EXTEND_GAP_A => println!("  Extend gap in A: {}", length),
            OP_EXTEND_GAP_B => println!("  Extend gap in B: {}", length),
            _ => println!("  Unknown operation: {} ({})", op, length),
        }
    }
*/

    // Reconstruct left aligned sequences
    right_edit_script.reverse();
    let (right_a, right_b) = reconstruct_alignment(&sequence_a, &sequence_b, right_edit_script, (seed_a_offset+1) as usize, (seed_b_offset+1) as usize);
    println!("Final Score = {}", total_score );
    println!("Final aligned range: A {}-{}", a_start, a_end );
    println!("Final aligned range: B {}-{}", b_start, b_end );

    // Combine left and right aligned segments
    println!("Sequence A: {}*{}", String::from_utf8(left_a).unwrap(), String::from_utf8(right_a).unwrap());
    println!("Sequence B: {}*{}", String::from_utf8(left_b).unwrap(), String::from_utf8(right_b).unwrap());
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_align_ex_left_right_extension() {
        // --- Step 1: Define sequences ---
        let sequence_a: Vec<u8> = b"TTCGTACGGCATCGAGTACC".to_vec();
        let sequence_b: Vec<u8> = b"TTCGTACGCAATCGAGTACC".to_vec();

        // Seed location
        let seed_a_offset = 15;
        let seed_b_offset = 15;

        // --- Step 2: Define scoring parameters ---
        let gap_open = 15;
        let gap_extend = 4;

        // Create a scoring matrix for DNA bases (A, C, G, T)
        let scoring_matrix = vec![
            //  A     C     G     T
            vec![9,   -15, -6,  -17], // A
            vec![-15, 10,  -15, -6 ], // C
            vec![-6,  -15, 10,  -15], // G
            vec![-17, -6,  -15, 9  ], // T
        ];

        // --- Step 3: Create the alignment structure ---
        let mut gap_align = BlastGapAlignStruct {
            dp_mem: Vec::new(),
            dp_mem_alloc: 0,
        };

        let scoring_params = BlastScoringParameters {
            scoring_matrix: scoring_matrix,
            gap_open,
            gap_extend,
            gap_x_dropoff: 50, // X-dropoff threshold
        };

        let mut total_score = 0;

        // --- Step 4: Run the left alignment extension ---
        let mut a_offset = 0;
        let mut b_offset = 0;

        let (left_best_score, left_edit_script) = align_ex(
            &sequence_a,
            &sequence_b,
            seed_a_offset as i32 + 1,
            seed_b_offset as i32 + 1,
            &mut a_offset,
            &mut b_offset,
            &mut gap_align,
            &scoring_params,
            true, // Reverse sequence direction (left extension)
        );

        // Save left extension results

        let a_start = seed_a_offset - a_offset + 1;
        let b_start = seed_b_offset - b_offset + 1;
        total_score += left_best_score;

        // --- Step 5: Run the right alignment extension ---
        a_offset = 0;
        b_offset = 0;

        let (right_best_score, right_edit_script) = align_ex(
            &sequence_a[(seed_a_offset as usize)..],
            &sequence_b[(seed_b_offset as usize)..],
            (sequence_a.len() - (seed_a_offset as usize) - 1) as i32,
            (sequence_b.len() - (seed_b_offset as usize) - 1) as i32,
            &mut a_offset,
            &mut b_offset,
            &mut gap_align,
            &scoring_params,
            false, // Forward sequence direction (right extension)
        );

        // Save right extension results

        let a_end = seed_a_offset + a_offset;
        let b_end = seed_b_offset + b_offset;
        total_score += right_best_score;

        // --- Step 6: Assertions ---
        assert_eq!(left_best_score, 105);
        assert_eq!(right_best_score, 38);
        assert_eq!(total_score, 143);
        assert_eq!(a_start, 0);
        assert_eq!(a_end, 19);
        assert_eq!(b_start, 0);
        assert_eq!(b_end, 19);
       
        let expected_left_edits = vec![
            OP_SUBSTITUTION,
            OP_SUBSTITUTION,
            OP_SUBSTITUTION,
            OP_SUBSTITUTION,
            OP_SUBSTITUTION,
            OP_SUBSTITUTION,
            OP_SUBSTITUTION,
            OP_SUBSTITUTION,
            OP_GAP_IN_B,
            OP_SUBSTITUTION,
            OP_SUBSTITUTION,
            OP_GAP_IN_A,
            OP_SUBSTITUTION,
            OP_SUBSTITUTION,
            OP_SUBSTITUTION,
            OP_SUBSTITUTION,
            OP_SUBSTITUTION,
        ];

        let expected_right_edits = vec![
            OP_SUBSTITUTION,
            OP_SUBSTITUTION,
            OP_SUBSTITUTION,
            OP_SUBSTITUTION,
        ];

        assert_eq!(left_edit_script, expected_left_edits);
        assert_eq!(right_edit_script, expected_right_edits);
    }

    #[test]
    fn test_align_ex_extended_gap() {
        // --- Step 1: Define sequences ---
        //                                          v
        let sequence_a: Vec<u8> = b"TGGATGGGGGAGGGAGTCCG".to_vec();
        let sequence_b: Vec<u8> = b"TGGATGGGGGTAGGAGGGAGTCCG".to_vec();
        //                                              ^

        // Seed location
        let seed_a_offset = 11;
        let seed_b_offset = 15;

        // --- Step 2: Define scoring parameters ---
        let gap_open = 20;
        let gap_extend = 5;

        // Create a scoring matrix for DNA bases (A, C, G, T)
        let scoring_matrix = vec![
            //  A     C     G     T
            vec![9,   -15, -6,  -17], // A
            vec![-15, 10,  -15, -6 ], // C
            vec![-6,  -15, 10,  -15], // G
            vec![-17, -6,  -15, 9  ], // T
        ];

        // --- Step 3: Create the alignment structure ---
        let mut gap_align = BlastGapAlignStruct {
            dp_mem: Vec::new(),
            dp_mem_alloc: 0,
        };

        let scoring_params = BlastScoringParameters {
            scoring_matrix: scoring_matrix,
            gap_open,
            gap_extend,
            gap_x_dropoff: 50, // X-dropoff threshold
        };

        let mut total_score = 0;

        // --- Step 4: Run the left alignment extension ---
        let mut a_offset = 0;
        let mut b_offset = 0;

        let (left_best_score, left_edit_script) = align_ex(
            &sequence_a,
            &sequence_b,
            seed_a_offset as i32 + 1,
            seed_b_offset as i32 + 1,
            &mut a_offset,
            &mut b_offset,
            &mut gap_align,
            &scoring_params,
            true, // Reverse sequence direction (left extension)
        );

        let a_start = seed_a_offset - a_offset + 1;
        let b_start = seed_b_offset - b_offset + 1;
        total_score += left_best_score;

        // --- Step 5: Run the right alignment extension ---
        a_offset = 0;
        b_offset = 0;

        let (right_best_score, right_edit_script) = align_ex(
            &sequence_a[(seed_a_offset as usize)..],
            &sequence_b[(seed_b_offset as usize)..],
            (sequence_a.len() - (seed_a_offset as usize) - 1) as i32,
            (sequence_b.len() - (seed_b_offset as usize) - 1) as i32,
            &mut a_offset,
            &mut b_offset,
            &mut gap_align,
            &scoring_params,
            false, // Forward sequence direction (right extension)
        );


        let a_end = seed_a_offset + a_offset;
        let b_end = seed_b_offset + b_offset;
        total_score += right_best_score;

        // --- Step 6: Assertions ---
        // 
        // 154 0.00 20.00 0.00 t5 1 20 (0) t6 1 24 (0)
        //  t5                     1 TGGATGGGGG----AGGGAGTCCG 20
        //                                     ----
        //  t6                     1 TGGATGGGGGTAGGAGGGAGTCCG 24

        assert_eq!(left_best_score, 76);
        assert_eq!(right_best_score, 78);
        assert_eq!(total_score, 154);
        assert_eq!(a_start, 0);
        assert_eq!(a_end, 19);
        assert_eq!(b_start, 0);
        assert_eq!(b_end, 23);

        let expected_left_edits = vec![
            OP_SUBSTITUTION,
            OP_SUBSTITUTION,
            OP_SUBSTITUTION,
            OP_SUBSTITUTION,
            OP_SUBSTITUTION,
            OP_SUBSTITUTION,
            OP_SUBSTITUTION,
            OP_SUBSTITUTION,
            OP_SUBSTITUTION,
            OP_SUBSTITUTION,
            OP_GAP_IN_A,
            OP_GAP_IN_A,
            OP_GAP_IN_A,
            OP_GAP_IN_A,
            OP_SUBSTITUTION,
            OP_SUBSTITUTION,
        ];

        let expected_right_edits = vec![
            OP_SUBSTITUTION,
            OP_SUBSTITUTION,
            OP_SUBSTITUTION,
            OP_SUBSTITUTION,
            OP_SUBSTITUTION,
            OP_SUBSTITUTION,
            OP_SUBSTITUTION,
            OP_SUBSTITUTION,
        ];

        assert_eq!(left_edit_script, expected_left_edits);
        assert_eq!(right_edit_script, expected_right_edits);
    }

    #[test]
    fn test_align_ex_alu1() {
        // --- Step 1: Define sequences ---

    // aluy
    let sequence_a: Vec<u8> = b"GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACGAGGTCAGGAGATCGAGACCATCCTGGCTAACACGGTGAAACCCCGTCTCTACTAAAAATACAAAAAATTAGCCGGGCGTGGTGGCGGGCGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGCAGGAGAATGGCGTGAACCCGGGAGGCGGAGCTTGCAGTGAGCCGAGATCGCGCCACTGCACTCCAGCCTGGGCGACAGAGCGAGACTCCGTCTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".to_vec();
    // example instance of an alu
    let sequence_b: Vec<u8> = b"GGCCGGGTGTGGTGGCTCATGCCTGCAATCCCAGCACTTTCAGGGGCCGAGGCGGGAGGATCGCTTGTGCTTAGGAGTTTGAGACTAGCCTGGGCAACATAGCAAGTCCCCGTCTCTACTTAAAAAAGAAAAAAAAATGAGCCAGATGTGGTGGCACGCAGCTGTAGTCCCAGCTACATGGAAGACTGAGGTGGTAGGATCACTTGAGCCAGGAGGTTGAGGCCACAGTGAGCCATGGTCACACCACTGCACTCCAGCCTGGGTGACAGAGTGAGACCCTGTCTCAAA".to_vec();
    
    // Seed location
    let seed_a_offset = 242;
    let seed_b_offset = 246; 

    // --- Step 2: Define scoring parameters ---
    let gap_open = 20;
    let gap_extend = 5;

    // Create a 5x5 scoring matrix for DNA bases (A, C, G, T, N)
    let scoring_matrix = vec![
      //   A   C   G   T   N
      vec![9, -15, -6, -17, -1], // A
      vec![-15, 10, -15, -6, -1], // C
      vec![-6, -15, 10, -15, -1], // G
      vec![-17, -6, -15, 9, -1], // T
      vec![-1, -1, -1, -1, -1], // N 
    ];


        // --- Step 3: Create the alignment structure ---
        let mut gap_align = BlastGapAlignStruct {
            dp_mem: Vec::new(),
            dp_mem_alloc: 0,
        };

        let scoring_params = BlastScoringParameters {
            scoring_matrix: scoring_matrix,
            gap_open,
            gap_extend,
            gap_x_dropoff: 200, // X-dropoff threshold
        };

        let mut total_score = 0;

        // --- Step 4: Run the left alignment extension ---
        let mut a_offset = 0;
        let mut b_offset = 0;

        let (left_best_score, left_edit_script) = align_ex(
            &sequence_a,
            &sequence_b,
            seed_a_offset as i32 + 1,
            seed_b_offset as i32 + 1,
            &mut a_offset,
            &mut b_offset,
            &mut gap_align,
            &scoring_params,
            true, // Reverse sequence direction (left extension)
        );

        let a_start = seed_a_offset - a_offset + 1;
        let b_start = seed_b_offset - b_offset + 1;

        // Reconstruct left aligned sequences
        let (left_a, left_b) = reconstruct_alignment(&sequence_a, &sequence_b, left_edit_script, a_start as usize, b_start as usize);

        total_score += left_best_score;

        // --- Step 5: Run the right alignment extension ---
        a_offset = 0;
        b_offset = 0;

        let (right_best_score, mut right_edit_script) = align_ex(
            &sequence_a[(seed_a_offset as usize)..],
            &sequence_b[(seed_b_offset as usize)..],
            (sequence_a.len() - (seed_a_offset as usize) - 1) as i32,
            (sequence_b.len() - (seed_b_offset as usize) - 1) as i32,
            &mut a_offset,
            &mut b_offset,
            &mut gap_align,
            &scoring_params,
            false, // Forward sequence direction (right extension)
        );
        
        let a_end = seed_a_offset + a_offset;
        let b_end = seed_b_offset + b_offset;

    // Reconstruct right aligned sequences
    right_edit_script.reverse();
    let (right_a, right_b) = reconstruct_alignment(
        &sequence_a,
        &sequence_b,
        right_edit_script,
        (seed_a_offset+1) as usize,  // Local start index for sliced sequence_a
        (seed_b_offset+1) as usize,  // Local start index for sliced sequence_b
    );

        total_score += right_best_score;

        // --- Step 6: Assertions ---
        // 
//1456 22.18 1.76 0.35 aluy 1 284 (27) aluex 1 288 (0)
//
//  aluy                   1 GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGA 50
//                                  i i         i     i              vi i      
//  aluex                  1 GGCCGGGTGTGGTGGCTCATGCCTGCAATCCCAGCACTTTCAGGGGCCGA 50
//
//  aluy                  51 GGCGGGCGGATCAC--GAGGTCAGGAGATCGAGACCATCCTGGCTAACAC 98
//                                 v     i -- v v i     v i     i v     vi    i
//  aluex                 51 GGCGGGAGGATCGCTTGTGCTTAGGAGTTTGAGACTAGCCTGGGCAACAT 100
//
//  aluy                  99 GGTGAAACCCCGTCTCTACT---AAAAATACAAAAAATTAGCCGGGCGTG 145
//                           i ii iv             ---    iv v       v    i ii   
//  aluex                101 AGCAAGTCCCCGTCTCTACTTAAAAAAGAAAAAAAAATGAGCCAGATGTG 150
//
//  aluy                 146 GTGGCGGGCGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGCAGGAGAAT 195
//                                iv  iv                vi  i  i      ii v  i  
//  aluex                151 GTGGCACGCAGCTGTAGTCCCAGCTACATGGAAGACTGAGGTGGTAGGAT 200
//
//  aluy                 196 GGCGTGAACCCGGGAGGCGGAGCTTGCAGTGAGCCGAGATCGCGCCACTG 245
//                           vi v   i  -i     iv   viii         iv i  i i      
//  aluex                201 CACTTGAGCC-AGGAGGTTGAGGCCACAGTGAGCCATGGTCACACCACTG 249
//
//  aluy                 246 CACTCCAGCCTGGGCGACAGAGCGAGACTCCGTCTCAAA 284
//                                         i       i     i i        
//  aluex                250 CACTCCAGCCTGGGTGACAGAGTGAGACCCTGTCTCAAA 288

// 1456	22.18	1.76	0.35	aluy	1	284	311	plus	aluex	1	288	288	27.27	26.77	20	6	43	GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCAC--GAGGTCAGGAGATCGAGACCATCCTGGCTAACACGGTGAAACCCCGTCTCTACT---AAAAATACAAAAAATTAGCCGGGCGTGGTGGCGGGCGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGCAGGAGAATGGCGTGAACCCGGGAGGCGGAGCTTGCAGTGAGCCGAGATCGCGCCACTGCACTCCAGCCTGGGCGACAGAGCGAGACTCCGTCTCAAA	GGCCGGGTGTGGTGGCTCATGCCTGCAATCCCAGCACTTTCAGGGGCCGAGGCGGGAGGATCGCTTGTGCTTAGGAGTTTGAGACTAGCCTGGGCAACATAGCAAGTCCCCGTCTCTACTTAAAAAAGAAAAAAAAATGAGCCAGATGTGGTGGCACGCAGCTGTAGTCCCAGCTACATGGAAGACTGAGGTGGTAGGATCACTTGAGCC-AGGAGGTTGAGGCCACAGTGAGCCATGGTCACACCACTGCACTCCAGCCTGGGTGACAGAGTGAGACCCTGTCTCAAA

    //println!("Sequence A: {}*{}", String::from_utf8(left_a).unwrap(), String::from_utf8(right_a).unwrap());
    //println!("Sequence B: {}*{}", String::from_utf8(left_b).unwrap(), String::from_utf8(right_b).unwrap());
        assert_eq!(left_best_score, 1125);
        assert_eq!(right_best_score, 331);
        assert_eq!(total_score, 1456);
        assert_eq!(a_start, 0);
        assert_eq!(a_end, 283);
        assert_eq!(b_start, 0);
        assert_eq!(b_end, 287);

        let a_aligned = String::from_utf8(left_a).unwrap() + &String::from_utf8(right_a).unwrap();
        let b_aligned = String::from_utf8(left_b).unwrap() + &String::from_utf8(right_b).unwrap();

        assert_eq!(a_aligned, "GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCAC--GAGGTCAGGAGATCGAGACCATCCTGGCTAACACGGTGAAACCCCGTCTCTACT---AAAAATACAAAAAATTAGCCGGGCGTGGTGGCGGGCGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGCAGGAGAATGGCGTGAACCCGGGAGGCGGAGCTTGCAGTGAGCCGAGATCGCGCCACTGCACTCCAGCCTGGGCGACAGAGCGAGACTCCGTCTCAAA");
        assert_eq!(b_aligned, "GGCCGGGTGTGGTGGCTCATGCCTGCAATCCCAGCACTTTCAGGGGCCGAGGCGGGAGGATCGCTTGTGCTTAGGAGTTTGAGACTAGCCTGGGCAACATAGCAAGTCCCCGTCTCTACTTAAAAAAGAAAAAAAAATGAGCCAGATGTGGTGGCACGCAGCTGTAGTCCCAGCTACATGGAAGACTGAGGTGGTAGGATCACTTGAGCC-AGGAGGTTGAGGCCACAGTGAGCCATGGTCACACCACTGCACTCCAGCCTGGGTGACAGAGTGAGACCCTGTCTCAAA");
    }
    

    #[test]
    fn test_align_ex_alu3_sub() {
        // --- Step 1: Define sequences ---
        //
    // aluy
    let sequence_a: Vec<u8> = b"GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACGAGGTCAGGAGATCGAGACCATCCTGGCTAACACGGTGAAACCCCGTCTCTACTAAAAATACAAAAAATTAGCCGGGCGTGGTGGCGGGCGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGCAGGAGAATGGCGTGAACCCGGGAGGCGGAGCTTGCAGTGAGCCGAGATCGCGCCACTGCACTCCAGCCTGGGCGACAGAGCGAGACTCCGTCTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".to_vec();
    // example instance of an alu
   let sequence_b: Vec<u8> = b"CCCCGTCTCTACTTAAAAAAGAAAAAAAAATGAGCCAGATGTG".to_vec();
    
    // Seed location
    let seed_a_offset = 108;
    let seed_b_offset = 3; 

    // --- Step 2: Define scoring parameters ---
    let gap_open = 20;
    let gap_extend = 5;

    // Create a 5x5 scoring matrix for DNA bases (A, C, G, T, N)
    let scoring_matrix = vec![
      //   A   C   G   T   N
      vec![9, -15, -6, -17, -1], // A
      vec![-15, 10, -15, -6, -1], // C
      vec![-6, -15, 10, -15, -1], // G
      vec![-17, -6, -15, 9, -1], // T
      vec![-1, -1, -1, -1, -1], // N 
    ];

    let x_dropoff = 200;
    let ( aligned_a, aligned_b, left_score, right_score, (a_start, a_end), (b_start, b_end) ) = align_sequences_with_seed(
        &sequence_a,
        &sequence_b,
        seed_a_offset,
        seed_b_offset,
        gap_open,
        gap_extend,
        scoring_matrix,
        x_dropoff,
    );

        // --- Step 6: Assertions ---
        // 
        //
// 205 17.50 7.50 0.00 aluy 106 145 (166) aluex2 1 43 (0)
//
//  aluy                 106 CCCCGTCTCTACT---AAAAATACAAAAAATTAGCCGGGCGTG 145
//                                        ---    iv v       v    i ii   
//  aluex2                 1 CCCCGTCTCTACTTAAAAAAGAAAAAAAAATGAGCCAGATGTG 43
//
//205	17.50	7.50	0.00	aluy	106	145	311	plus	aluex2	1	43	43	20.14	20.14	3	1	4	CCCCGTCTCTACT---AAAAATACAAAAAATTAGCCGGGCGTG	CCCCGTCTCTACTTAAAAAAGAAAAAAAAATGAGCCAGATGTG
 
        assert_eq!(left_score, 40);
        assert_eq!(right_score, 165);
        assert_eq!(a_start, 105);
        assert_eq!(a_end, 144);
        assert_eq!(b_start, 0);
        assert_eq!(b_end, 42);

        assert_eq!(aligned_a, "CCCCGTCTCTACT---AAAAATACAAAAAATTAGCCGGGCGTG", "sequence_a alignment string mismatch");
        assert_eq!(aligned_b, "CCCCGTCTCTACTTAAAAAAGAAAAAAAAATGAGCCAGATGTG", "sequence_b alignment string mismatch");
    }
    #[test]
    fn test_align_ex_complex_left_aligned() {
        // --- Step 1: Define sequences ---
        //
     //// aluy
    let sequence_a: Vec<u8> = b"GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACGAGGTCAGGAGATCGAGACCATCCTGGCTAACACGGTGAAACCCCGTCTCTACTAAAAATACAAAAAATTAGCCGGGCGTGGTGGCGGGCGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGCAGGAGAATGGCGTGAACCCGGGAGGCGGAGCTTGCAGTGAGCCGAGATCGCGCCACTGCACTCCAGCCTGGGCGACAGAGCGAGACTCCGTCTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".to_vec();
    // example instance of an alu
    // 220	12.50	10.00	0.00	aluy	106	145	311	plus	aluex2	1	44	44	13.79	13.79	2	2	3	CCCCGTCTCTACT---AAAAATACAAAA-AATTAGCCGGGCGTG	CCGCGTCTCTACTTAAAAAAATAAAAAAAAATTAGCCAGATGTG
    let sequence_b: Vec<u8> = b"CCGCGTCTCTACTTAAAAAAATAAAAAAAAATTAGCCAGATGTG".to_vec();
    
    // Seed location
    let seed_a_offset = 129;
    let seed_b_offset = 28; 

    // --- Step 2: Define scoring parameters ---
    let gap_open = 20;
    let gap_extend = 5;

    // Create a 5x5 scoring matrix for DNA bases (A, C, G, T, N)
    let scoring_matrix = vec![
      //   A   C   G   T   N
      vec![9, -15, -6, -17, -1], // A
      vec![-15, 10, -15, -6, -1], // C
      vec![-6, -15, 10, -15, -1], // G
      vec![-17, -6, -15, 9, -1], // T
      vec![-1, -1, -1, -1, -1], // N 
    ];

    let x_dropoff = 200;
    let ( aligned_a, aligned_b, left_score, right_score, (a_start, a_end), (b_start, b_end) ) = align_sequences_with_seed(
        &sequence_a,
        &sequence_b,
        seed_a_offset,
        seed_b_offset,
        gap_open,
        gap_extend,
        scoring_matrix,
        x_dropoff,
    );

        // --- Step 6: Assertions ---
        // 
        //
//  aluy                 106 CCCCGTCTCTACT---AAAAATACAAAA-AATTAGCCGGGCGTG 145
//                             v          ---       v    -        i ii   
//  aluex2                 1 CCGCGTCTCTACTTAAAAAAATAAAAAAAAATTAGCCAGATGTG 44

//220 12.50 10.00 0.00  aluy  106 145 311 plus  aluex2  1 44  44  13.79 13.79 2 2 3 CCCCGTCTCTACT---AAAAATACAAAA-AATTAGCCGGGCGTG  CCGCGTCTCTACTTAAAAAAATAAAAAAAAATTAGCCAGATGTG       
 
        assert_eq!(left_score, 124);
        assert_eq!(right_score, 96);
        assert_eq!(a_start, 105);
        assert_eq!(a_end, 144);
        assert_eq!(b_start, 0);
        assert_eq!(b_end, 43);

        assert_eq!(aligned_a, "CCCCGTCTCTACT---AAAAATACAAAA-AATTAGCCGGGCGTG", "sequence_a alignment string mismatch");
        assert_eq!(aligned_b, "CCGCGTCTCTACTTAAAAAAATAAAAAAAAATTAGCCAGATGTG", "sequence_b alignment string mismatch");
    }

}


