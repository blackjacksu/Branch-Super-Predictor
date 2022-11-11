// File: Two level branch predictor
// Date:
// Institute: ECE, Texas A&M University
// Author: TING-WEI SU (willytwsu@tamu.edu)

// Macro Pattern History Table (PHT)
#define J_BITS 6    // Use the last j-bits in Branch address to index the PHT
#define PATT_HIST_TBL_NUM   (1 << J_BITS)
#define PATT_HIST_TBL_NUM_MAX ((1 << J_BITS) - 1)
#define J_BITS_MASK     PATT_HIST_TBL_NUM_MAX

// Macro Branch History Shift Register (BHSR)
#define I_BITS  0
#define K_BITS  8    // Use last k-bits in BHSR to index the PHT
#define PATT_HIST_TBL_ENTRY_NUM     (1 << K_BITS)
#define PATT_HIST_TBL_ENTRY_NUM_MAX ((1 << K_BITS) - 1)
#define K_BITS_MASK     PATT_HIST_TBL_ENTRY_NUM_MAX

enum pred {
    strong_nottaken,
    weak_nottaken,
    weak_taken,
    strong_taken,
};

class pattern_history_table {
    private:
        enum pred two_bit_ctr[PATT_HIST_TBL_NUM][PATT_HIST_TBL_ENTRY_NUM];

    public:
        pattern_history_table() {}
        
        ~pattern_history_table() {}

        bool get_prediction(unsigned long long pht_idx, unsigned long long entry_idx)
        {
            return two_bit_ctr[pht_idx][entry_idx] >= 2;
        }

        void update_counter_taken(unsigned long long pht_idx, unsigned long long entry_idx)
        {
            switch (two_bit_ctr[pht_idx][entry_idx])
            {
                case strong_nottaken: // Strongly not taken
                    two_bit_ctr[pht_idx][entry_idx] = weak_nottaken;
                    break;
                case weak_nottaken: // Weakly not taken
                    two_bit_ctr[pht_idx][entry_idx] = weak_taken;
                    break;
                case weak_taken: // Weakly taken
                    two_bit_ctr[pht_idx][entry_idx] = strong_taken;
                    break;
                case strong_taken: // Strongly taken
                    two_bit_ctr[pht_idx][entry_idx] = strong_taken;
                    break;
                default:
                    // The first time update taken
                    two_bit_ctr[pht_idx][entry_idx] = weak_taken;
                    break;
            }
        }

        void update_counter_nottaken(unsigned long long pht_idx, unsigned long long entry_idx)
        {
            switch (two_bit_ctr[pht_idx][entry_idx])
            {
                case strong_nottaken: // Strongly not taken
                    two_bit_ctr[pht_idx][entry_idx] = strong_nottaken;
                    break;
                case weak_nottaken: // Weakly not taken
                    two_bit_ctr[pht_idx][entry_idx] = strong_nottaken;
                    break;
                case weak_taken: // Weakly taken
                    two_bit_ctr[pht_idx][entry_idx] = weak_nottaken;
                    break;
                case strong_taken: // Strongly taken
                    two_bit_ctr[pht_idx][entry_idx] = weak_taken;
                    break;
                default: 
                    // The first time update not taken
                    two_bit_ctr[pht_idx][entry_idx] = weak_nottaken;
                    break;
            }
        }
};

class twolevel_predictor : public branch_predictor
{
    private:
        class pattern_history_table pht_tbl;
        unsigned long long br_hist_shift_reg;
        unsigned int bhsr_k_idx;
        unsigned int br_addr_j_idx;

    public:
    	branch_update u;
	    branch_info bi;

        twolevel_predictor()
        {
            br_hist_shift_reg = 0;
        }

        ~twolevel_predictor() {}

        // PREDICTION
	    branch_update *predict (branch_info & b)
	    {
            bi = b;
            if (b.br_flags & BR_CONDITIONAL) 
            {
                br_addr_j_idx = bi.address & J_BITS_MASK;
                bhsr_k_idx = br_hist_shift_reg & K_BITS_MASK;
			    u.direction_prediction (pht_tbl.get_prediction(br_addr_j_idx, bhsr_k_idx));
		    } 
            else 
            {
			    u.direction_prediction (true);
		    }
		    u.target_prediction (0);
		    return &u;
        }

        void update (branch_update * u, bool taken, unsigned int target)
	    {
            if (bi.br_flags & BR_CONDITIONAL) 
            {
                br_addr_j_idx = bi.address & J_BITS_MASK;
                br_hist_shift_reg = (br_hist_shift_reg << 1) | taken;
                bhsr_k_idx = br_hist_shift_reg & K_BITS_MASK;

                if (taken)
                {
                    pht_tbl.update_counter_taken(br_addr_j_idx, bhsr_k_idx);
                }
                else
                {
                    pht_tbl.update_counter_nottaken(br_addr_j_idx, bhsr_k_idx);
                }
		    }
        }
};