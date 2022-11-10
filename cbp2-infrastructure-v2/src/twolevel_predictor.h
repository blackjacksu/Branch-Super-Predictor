// File: Two level branch predictor
// Date:
// Institute: ECE, Texas A&M University
// Author: TING-WEI SU (willytwsu@tamu.edu)

#define BR_HISTORY_REG_LEN 16
#define BR_HISTORY_REG_MAX ((1 << BR_HISTORY_REG_LEN) - 1)
#define PATTERN_HISTORY_ENTRY (1 << BR_HISTORY_REG_LEN)


enum pred {
    strong_nottaken,
    weak_nottaken,
    weak_taken,
    strong_taken,
};

class pattern_history_table {
    private:
        unsigned long long num_PHT_entry;
        enum pred* two_bit_counter;


    public:
        pattern_history_table()
        {
            num_PHT_entry = PATTERN_HISTORY_ENTRY;
            two_bit_counter = new enum pred[num_PHT_entry];
        }
        
        ~pattern_history_table()
        {
            delete two_bit_counter;
        }

        void update_counter_taken(unsigned long long index)
        {
            switch (two_bit_counter[index])
            {
                case strong_nottaken: // Strongly not taken
                    two_bit_counter[index] = weak_nottaken;
                    break;
                case weak_nottaken: // Weakly not taken
                    two_bit_counter[index] = weak_taken;
                    break;
                case weak_taken: // Weakly taken
                    two_bit_counter[index] = strong_taken;
                    break;
                case strong_taken: // Strongly taken
                    two_bit_counter[index] = strong_taken;
                    break;
                default:
                    two_bit_counter[index] = weak_taken;
                    break;
            }
        }

        void update_counter_nottaken(unsigned long long index)
        {
            switch (two_bit_counter[index])
            {
                case strong_nottaken: // Strongly not taken
                    two_bit_counter[index] = strong_nottaken;
                    break;
                case weak_nottaken: // Weakly not taken
                    two_bit_counter[index] = strong_nottaken;
                    break;
                case weak_taken: // Weakly taken
                    two_bit_counter[index] = weak_nottaken;
                    break;
                case strong_taken: // Strongly taken
                    two_bit_counter[index] = weak_taken;
                    break;
                default:
                    two_bit_counter[index] = weak_nottaken;
                    break;
            }
        }

        bool get_prediction(unsigned long long index)
        {
            return two_bit_counter[index] > 2;
        }
};

class twolevel_predictor : public branch_predictor
{
    private:
        class pattern_history_table pattern_history_table;
        unsigned long long branch_history_shift_register;

    public:
    	branch_update u;
	    branch_info bi;

        twolevel_predictor()
        {
            branch_history_shift_register = 0;
        }

        ~twolevel_predictor()
        {

        }

        // PREDICTION
	    branch_update *predict (branch_info & b)
	    {
            bi = b;
            if (b.br_flags & BR_CONDITIONAL) {
			    // u.index = 
				//   (history << (TABLE_BITS - HISTORY_LENGTH)) 
				// ^ (b.address & ((1<<TABLE_BITS)-1));
			    u.direction_prediction (pattern_history_table.get_prediction(branch_history_shift_register));
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
                branch_history_shift_register = ((branch_history_shift_register << 1) | taken) & BR_HISTORY_REG_MAX; // Use BR_HISTORY_REG_MAX to Mask the length of bhsr to prevent seg fault in pht
                if (taken)
                {
                    pattern_history_table.update_counter_taken(branch_history_shift_register);
                }
                else
                {
                    pattern_history_table.update_counter_nottaken   (branch_history_shift_register);
                }
		    }
        }
};