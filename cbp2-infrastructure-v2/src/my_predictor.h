/* 
Code has been largely inspired	by the tagged PPM predictor simulator from Pierre Michaud, the OGEHL predictor simulator from by André Seznec, the TAGE predictor simulator from	André Seznec and Pierre Michaud



*/
// my_predictor.h
#include <inttypes.h>
#include <math.h>


//a limit predictor for 256Kbits: no AHEAD pipelining, 13 components
#define NHIST 18		//12 tagged components                               // 要改
#define LOGB 14			// log of the number of entries in the base bimodal predictor
#define HYSTSHIFT 2		// sharing an hysteris bit between 4 bimodal predictor entries
#define LOGG (LOGB-3)		// base 2 logarithm of number of entries	on each tagged component
#define TBITS 7			// minimum tag width (shortest history length table)
#define MINHIST 4		// shortest history length
#define MAXHIST 1280		// longest history lngth

#define LOGL 8			//256 entries loop predictor
#define WIDTHNBITERLOOP 14
//Total storage for the submitted predictor
//TAGE
//Table T0: 20Kbits (16K prediction + 4K hysteresis)
//Tables T1 and T2: 12Kbits each; Tables T3 and T4: 26 Kbits each; Table T5: 28Kbits; Table	T6: 30Kbits; Table T7: 16 Kbits	;	Tables T8 and	T9: 17 Kbits each, Table T10: 18 Kbits;	Table T11: 9,5 Kbits; Table T12 10 Kbits
// Total of storage for the TAGE tables: 239,5 Kbits
//Loop predictor: 256 entries * 52 bits= 13 Kbits
// Total of storage for the TAGE tables: 241,5 Kbits + 13 Kbits = 260,608 bits
//Extra storage: 2*640 history bits + 2*16 path history bits + 4 bits for	USE_ALT_ON_NA + 19 bits for the TICK counter + 2 bits for Seed in the "pseudo-random" counter= 1,337 bits + 7bits on the WITHLOOP counter = 1,344 bits
//Total storage= 260,608 + 1,344 = 261,952 bits



#define CWIDTH 3		// predictor counter width on the tagged tables
typedef uint32_t address_t;

#define BUFFERHIST 128 * 1024
//size of the history buffer: to allow fast management we use a very large buffer //we replace a global shift by a pointer management




// utility class for index computation

// this is the cyclic shift register for folding 
// a long global history into a smaller number of bits; see P. Michaud's PPM-like predictor at CBP-1
class folded_history
{
public:
	unsigned comp;
	int CLENGTH;
	int OLENGTH;
	int OUTPOINT;

	folded_history ()
	{
	}

	void init (int original_length, int compressed_length)
	{
		comp = 0;
		OLENGTH = original_length;
		CLENGTH = compressed_length;
		OUTPOINT = OLENGTH % CLENGTH;
	}

	void update (uint8_t * h)
	{
		comp = (comp << 1) | h[0];
		comp ^= h[OLENGTH] << OUTPOINT;
		comp ^= (comp >> CLENGTH);
		comp &= (1 << CLENGTH) - 1;
	}
};


class TAGE : public branch_predictor
{
public:
	branch_update u;
	branch_info bi;

	//	class lentry			//our predictor entry

	class bentry			// TAGE bimodal table entry	
	{										 // base predictor
	public:
		int8_t hyst;
		int8_t pred;
		bentry ()
		{
			pred = 0;
			hyst = 1;
		}
	};

	class gentry			// TAGE global table entry
	{										 // other banks
	public:
		int8_t ctr;
		uint16_t tag;
		int8_t u;
		gentry ()
		{

			ctr = 0;
			tag = 0;
			u = 0;
		}
	};
	
	// "Use alternate prediction on newly allocated":	a 4-bit counter	to determine whether the	newly allocated entries should be considered as	valid or not for delivering	the prediction
	int USE_ALT_ON_NA;			
	//control counter for the smooth resetting of useful counters
	int TICK, LOGTICK;		

	// use a path history as on	the OGEHL predictor
	int phist;			

	//for managing global history	
	uint8_t *GHIST;
	uint8_t *ghist;
	int ptghist;
	// path history including kernel activity
	int phistos;				
	// for managing global history including kernel activity
	uint8_t *GHISTOS;
	uint8_t *ghistos;
	int ptghistos;	
	folded_history ch_i[NHIST + 1];	//utility for computing TAGE indices
	folded_history ch_t[2][NHIST + 1];	//utility for computing TAGE tags	
	folded_history ch_ios[NHIST + 1];	//utility for computing TAGE indices on kernel branches
	folded_history ch_tos[2][NHIST + 1];	//utility for computing TAGE tags on kernel branches
	//lentry *ltable;		//loop predictor table
	bentry *btable;		//bimodal TAGE table
	gentry *gtable[NHIST + 1];	// tagged TAGE tables	
	int TB[NHIST + 1];		// tag width for the different tagged tables
	int m[NHIST + 1];		// used for storing the history lengths
	int logg[NHIST + 1];		// log of number entries of the different tagged tables
	int GI[NHIST + 1];		// indexes to the different tables are computed only once	
	int GTAG[NHIST + 1];		// tags for the different tables are computed only once	
	int BI;			// index of the bimodal table
	int Seed;			// for the pseudo-random number generator
	bool pred_taken;		// prediction
	bool alttaken;		// alternate	TAGEprediction
	bool tage_pred;		// TAGE prediction
	int HitBank;			// longest matching bank
	int AltBank;			// alternate matching bank	
	//bool predloop;		// our predictor parameter{}	
	TAGE (void)
	{
		USE_ALT_ON_NA = 0;
		Seed = 0;
		LOGTICK = 19;		//lg smooth resetting
		TICK = (1 << (LOGTICK - 1));	//initog of the period for useful taialize the resetting counter to the half of the period
		phist = 0;
		phistos = 0;
		GHIST = (uint8_t *) malloc(BUFFERHIST);
		GHISTOS = (uint8_t *) malloc(BUFFERHIST);
		ghist = GHIST;
		ghistos = GHISTOS;

		for (int i = 0; i < BUFFERHIST; i++)
			ghist[0] = 0;

		for (int i = 0; i < BUFFERHIST; i++)
			ghistos[0] = 0;

		ptghist = 0;
		ptghistos = 0;

		// computes the geometric history lengths	 
		m[1] = MINHIST;
		m[NHIST] = MAXHIST;
		for (int i = 2; i <= NHIST; i++)
		{
			m[i] = (int) (((double) MINHIST *
				pow ((double) (MAXHIST) / (double) MINHIST,
				(double) (i - 1) / (double) ((NHIST - 1)))) + 0.5);
		}

		//widths of the partial tags
		TB[1] = TBITS;
		TB[2] = TBITS;
		TB[3] = TBITS + 1;
		TB[4] = TBITS + 1;
		TB[5] = TBITS + 2;
		TB[6] = TBITS + 3;
		TB[7] = TBITS + 4;
		TB[8] = TBITS + 5;
		TB[9] = TBITS + 5;
		TB[10] = TBITS + 6;
		TB[11] = TBITS + 7;
		TB[12] = TBITS + 8;
    	TB[13] = TBITS + 8;
    	TB[14] = TBITS + 9;
    	TB[15] = TBITS + 9;
    	TB[16] = TBITS + 10;
    	TB[17] = TBITS + 10;
    	TB[18] = TBITS + 11;

		// log2 of number entries in the tagged components
		for (int i = 1; i <= 2; i++)
			logg[i] = LOGG - 1;
		for (int i = 3; i <= 9; i++)
			logg[i] = LOGG;
		for (int i = 10; i <= 12; i++)
			logg[i] = LOGG - 1;
		for (int i = 13; i <= 14; i++)
			logg[i] = LOGG - 1;
 		for (int i = 14; i <=16 ; i++)
			logg[i] = LOGG;
 		for (int i = 17; i <=18 ; i++)
			logg[i] = LOGG -2;

		//initialisation of index and tag computation functions
		for (int i = 1; i <= NHIST; i++)
		{
			ch_i[i].init (m[i], (logg[i]));
			ch_t[0][i].init (ch_i[i].OLENGTH, TB[i]);
			ch_t[1][i].init (ch_i[i].OLENGTH, TB[i] - 1);
		}
		for (int i = 1; i <= NHIST; i++)
		{
			ch_ios[i].init (m[i], (logg[i]));
			ch_tos[0][i].init (ch_i[i].OLENGTH, TB[i]);
			ch_tos[1][i].init (ch_i[i].OLENGTH, TB[i] - 1);
		}
		//allocation of the loop predictor table
		//		ltable = new lentry[1 << LOGL];
		//allocation of the predictor tables
		btable = new bentry[1 << LOGB];
		for (int i = 1; i <= NHIST; i++)
		{
			gtable[i] = new gentry[1 << (logg[i])];
		}
	}


		// index function for the bimodal table

	int bindex (address_t pc)
	{
		return ((pc) & ((1 << (LOGB)) - 1));
	}

	// index function for the	4-way associative loop predictor
	int lindex (address_t pc)
	{
		return ((pc & ((1 << (LOGL - 2)) - 1)) << 2);
	}

	// the index functions for the tagged tables uses path history as in the OGEHL predictor
	//F serves to mix path history
	int F (int A, int size, int bank)
	{
		int A1, A2;

		A = A & ((1 << size) - 1);
		A1 = (A & ((1 << logg[bank]) - 1));
		A2 = (A >> logg[bank]);
		A2 =
			((A2 << bank) & ((1 << logg[bank]) - 1)) + (A2 >> (logg[bank] - bank));
		A = A1 ^ A2;
		A = ((A << bank) & ((1 << logg[bank]) - 1)) + (A >> (logg[bank] - bank));
		return (A);
	}

	// gindex computes a full hash of pc, ghist and phist
	int gindex (address_t pc, int bank)
	{
		int index;
		int M = (m[bank] > 16) ? 16 : m[bank];
		index =
			pc ^ (pc >> (abs (logg[bank] - bank) + 1)) ^
			ch_i[bank].comp ^ F (phist, M, bank);

		return (index & ((1 << (logg[bank])) - 1));
	}

	//	tag computation
	uint16_t gtag (address_t pc, int bank)
	{
		int tag = pc ^ ch_t[0][bank].comp ^ (ch_t[1][bank].comp << 1);

		return (tag & ((1 << TB[bank]) - 1));
	}

	//index computation for kernel branchs 
	int gindexos (address_t pc, int bank)
	{
		int index;
		int M = (m[bank] > 16) ? 16 : m[bank];
		index =
			pc ^ (pc >> (abs (logg[bank] - bank) + 1)) ^ ch_ios[bank].
			comp ^ F (phistos, M, bank);

		return (index & ((1 << (logg[bank])) - 1));
	}

	//	tag computation for kernel branch
	uint16_t gtagos (address_t pc, int bank)
	{
		int tag = pc ^ ch_tos[0][bank].comp ^ (ch_tos[1][bank].comp << 1);
		return (tag & ((1 << TB[bank]) - 1));
	}

	// up-down saturating counter
	void ctrupdate (int8_t & ctr, bool taken, int nbits)
	{
		if (taken)
		{
			if (ctr < ((1 << (nbits - 1)) - 1))
				ctr++;
		}
		else
		{
			if (ctr > -(1 << (nbits - 1)))
				ctr--;
		}
	}
	bool getbim (address_t pc)
	{
		return (btable[BI].pred > 0);
	}

	// update	the bimodal predictor: a hysteresis bit is shared among 4 prediction bits
	void baseupdate (address_t pc, bool Taken)
	{
		int inter = (btable[BI].pred << 1) + btable[BI >> HYSTSHIFT].hyst;
		if (Taken)
		{
			if (inter < 3)
				inter += 1;
		}
		else if (inter > 0)
		{
			inter--;
		}

		btable[BI].pred = inter >> 1;
		btable[BI >> HYSTSHIFT].hyst = (inter & 1);
	}

	//our predictor and its update function 
 

	//just a simple pseudo random number generator: a 2-bit counter, used to avoid ping-pong phenomenom on tagged entry allocations

	int MYRANDOM ()
	{
		Seed++;
		return (Seed & 3);
	};

	// shifting the global history:	we manage the history in a big table in order to reduce simulation time
	void updateghist (uint8_t * &h, bool dir, uint8_t * tab, int &PT)
	{
		if (PT == 0)
		{
			for (int i = 0; i < MAXHIST; i++)
			{
				tab[BUFFERHIST - MAXHIST + i] = tab[i];
			}
			PT = BUFFERHIST - MAXHIST;
			h = &tab[PT];
		}
		PT--;
		h--;
		h[0] = (dir) ? 1 : 0;
	}

	// PREDICTION
	branch_update *predict (branch_info & b)
	{
		bi = b;

		int pc = b.address;

		bool KERNELMODE = ((pc & 0xc0000000) == 0xc0000000);

		if (b.br_flags & BR_CONDITIONAL)
		{
			// TAGE prediction

			// computes the table addresses and the partial tags
			if (!KERNELMODE)
			{
				for (int i = 1; i <= NHIST; i++)
				{
					GI[i] = gindex (pc, i);
					GTAG[i] = gtag (pc, i);
				}
			}
			else
			{
				for (int i = 1; i <= NHIST; i++)
				{
					GI[i] = gindexos (pc, i);
					GTAG[i] = gtagos (pc, i);
				}
			}
	
			BI = pc & ((1 << LOGB) - 1);

			HitBank = 0;
			AltBank = 0;
			
			//Look for the bank with longest matching history
			for (int i = NHIST; i > 0; i--)
			{
				if (gtable[i][GI[i]].tag == GTAG[i])
				{
					HitBank = i;
					break;
				}
			}
			//Look for the alternate bank
			for (int i = HitBank - 1; i > 0; i--)
			{
				if (gtable[i][GI[i]].tag == GTAG[i])
					if ((USE_ALT_ON_NA < 0)
						|| (abs (2 * gtable[i][GI[i]].ctr + 1) > 1))
					{
						AltBank = i;
						break;
					}
			}
			//computes the prediction and the alternate prediction
			if (HitBank > 0)
			{
				if (AltBank > 0)
					alttaken = (gtable[AltBank][GI[AltBank]].ctr >= 0);
				else
					alttaken = getbim (pc);
			//if the entry is recognized as a newly allocated entry and 
			//USE_ALT_ON_NA is positive	use the alternate prediction
				if ((USE_ALT_ON_NA < 0)
					|| (abs (2 * gtable[HitBank][GI[HitBank]].ctr + 1) > 1))
					tage_pred = (gtable[HitBank][GI[HitBank]].ctr >= 0);
				else
					tage_pred = alttaken;
			}
			else
			{
				alttaken = getbim (pc);
				tage_pred = alttaken;
			}
			//end TAGE prediction

			//pred_taken = ((WITHLOOP >= 0) && (LVALID)) ? predloop : tage_pred;
			pred_taken = tage_pred;

		}
		u.direction_prediction (pred_taken);
		u.target_prediction (0);
		return &u;
	}

	// PREDICTOR UPDATE
	void update (branch_update * u, bool taken, unsigned int target)
	{
		int NRAND = MYRANDOM ();
		address_t pc = bi.address;
		if (bi.br_flags & BR_CONDITIONAL)
		{
			// update our predictor here
			// TAGE UPDATE	
			// try to allocate a	new entries only if prediction was wrong
			bool ALLOC = ((tage_pred != taken) & (HitBank < NHIST));
			if (HitBank > 0)
			{
				// Manage the selection between longest matching and alternate matching
				// for "pseudo"-newly allocated longest matching entry
				bool LongestMatchPred = (gtable[HitBank][GI[HitBank]].ctr >= 0);
				bool PseudoNewAlloc =
					(abs (2 * gtable[HitBank][GI[HitBank]].ctr + 1) <= 1);
				// an entry is considered as newly allocated if its prediction counter is weak
				if (PseudoNewAlloc)
				{
					if (LongestMatchPred == taken)
						ALLOC = false;
					// if it was delivering the correct prediction, no need to allocate a new entry
					//even if the overall prediction was false

					if (LongestMatchPred != alttaken)
					{
						if (alttaken == taken)
						{
							if (USE_ALT_ON_NA < 7)
								USE_ALT_ON_NA++;
						}
						else if (USE_ALT_ON_NA > -8)
							USE_ALT_ON_NA--;
					}

					if (USE_ALT_ON_NA >= 0)
						tage_pred = LongestMatchPred;
				}
			}

			if (ALLOC)
			{
				// is there some "unuseful" entry to allocate
				int8_t min = 1;
				for (int i = NHIST; i > HitBank; i--)
					if (gtable[i][GI[i]].u < min)
						min = gtable[i][GI[i]].u;

				// we allocate an entry with a longer history
				//to	avoid ping-pong, we do not choose systematically the next entry, but among the 3	next	 entries
				int Y = NRAND & ((1 << (NHIST - HitBank - 1)) - 1);
				int X = HitBank + 1;
				if (Y & 1)
				{
					X++;
					if (Y & 2)
						X++;
				}
				//NO ENTRY AVAILABLE:	ENFORCES ONE TO BE AVAILABLE 
				if (min > 0)
					gtable[X][GI[X]].u = 0;

				//Allocate only	one entry

				for (int i = X; i <= NHIST; i += 1)
				{
					if ((gtable[i][GI[i]].u == 0))
					{
						gtable[i][GI[i]].tag = GTAG[i];
						gtable[i][GI[i]].ctr = (taken) ? 0 : -1;
						gtable[i][GI[i]].u = 0;
						break;
					}
				}
			}
			//periodic reset of u: reset is not complete but bit by bit
			TICK++;
			if ((TICK & ((1 << LOGTICK) - 1)) == 0)
			{
				// reset least significant bit
				// most significant bit becomes least significant bit
				for (int i = 1; i <= NHIST; i++)
					for (int j = 0; j < (1 << logg[i]); j++)
						gtable[i][j].u = gtable[i][j].u >> 1;
			}

			if (HitBank > 0)
			{
				ctrupdate (gtable[HitBank][GI[HitBank]].ctr, taken, CWIDTH);
				//if the provider entry is not certified to be useful also update the alternate prediction
				if (gtable[HitBank][GI[HitBank]].u == 0)
				{
					if (AltBank > 0)
						ctrupdate (gtable[AltBank][GI[AltBank]].ctr, taken, CWIDTH);
					if (AltBank == 0)
						baseupdate (pc, taken);
				}
			}
			else
				baseupdate (pc, taken);

			// update the u counter
			if (tage_pred != alttaken)
			{
				if (tage_pred == taken)
				{
					if (gtable[HitBank][GI[HitBank]].u < 3)
						gtable[HitBank][GI[HitBank]].u++;
				}
				else
				{
					if (USE_ALT_ON_NA < 0)
						if (gtable[HitBank][GI[HitBank]].u > 0)
							gtable[HitBank][GI[HitBank]].u--;
				}
			}
		
			//END PREDICTOR UPDATE
		}

		//	UPDATE HISTORIES	

		//check user and kernel mode to detect transition:

		bool KERNELMODE = ((pc & 0xc0000000) == 0xc0000000);

		bool TAKEN = ((!(bi.br_flags & BR_CONDITIONAL)) | (taken));
		bool PATHBIT = (bi.address & 1);

		if (!KERNELMODE)
		{
			//update user history
			updateghist (ghist, TAKEN, GHIST, ptghist);
			phist = (phist << 1) + PATHBIT;

			phist = (phist & ((1 << 16) - 1));
			//prepare next index and tag computations for user branchs 
			for (int i = 1; i <= NHIST; i++)
			{
				ch_i[i].update (ghist);
				ch_t[0][i].update (ghist);
				ch_t[1][i].update (ghist);
			}
		}

		// always update kernel history
		updateghist (ghistos, TAKEN, GHISTOS, ptghistos);
		phistos = (phistos << 1) + PATHBIT;
		phistos = (phistos & ((1 << 16) - 1));
		//prepare next index and tag computations for	kernel branchs
		for (int i = 1; i <= NHIST; i++)
		{
			ch_ios[i].update (ghistos);
			ch_tos[0][i].update (ghistos);
			ch_tos[1][i].update (ghistos);
		}
		//END UPDATE HISTORIES
	}
};

////////////////////////////////////
// Two level branch predictor
////////////////////////////////////
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

#define STRONG_NOTTAKEN     0
#define MAJOR_NOTTAKEN      1
#define MINOR_NOTTAKEN      2
#define WEAK_NOTTAKEN       3
#define WEAK_TAKEN          4
#define MINOR_TAKEN         5
#define MAJOR_TAKEN         6
#define STRONG_TAKEN        7

enum pred {
    strong_nottaken,
    medium_nottaken,
    mediumrare_nottaken,
    weak_nottaken,
    weak_taken,
    mediumrare_taken,
    medium_taken,
    strong_taken,
};

class pattern_history_table {
    private:
        unsigned int two_bit_ctr[PATT_HIST_TBL_NUM][PATT_HIST_TBL_ENTRY_NUM];

    public:
        pattern_history_table() {}
        
        ~pattern_history_table() {}

        bool get_prediction(unsigned long long pht_idx, unsigned long long entry_idx)
        {
            return two_bit_ctr[pht_idx][entry_idx] >= weak_taken;
        }

        void update_counter_taken(unsigned long long pht_idx, unsigned long long entry_idx)
        {
            switch (two_bit_ctr[pht_idx][entry_idx])
            {
                case STRONG_NOTTAKEN: // Strongly not taken lower bound
                case MAJOR_NOTTAKEN:
                case MINOR_NOTTAKEN:
                case WEAK_NOTTAKEN: 
                case WEAK_TAKEN:   
                case MINOR_TAKEN:   
                case MAJOR_TAKEN:   
                    two_bit_ctr[pht_idx][entry_idx]++;
                    break;
                case STRONG_TAKEN: // Strongly taken upper bound
                    two_bit_ctr[pht_idx][entry_idx] = STRONG_TAKEN;
                    break;
                default:
                    // The first time update taken
                    two_bit_ctr[pht_idx][entry_idx] = WEAK_TAKEN;
                    break;
            }
        }

        void update_counter_nottaken(unsigned long long pht_idx, unsigned long long entry_idx)
        {
            switch (two_bit_ctr[pht_idx][entry_idx])
            {
                case STRONG_NOTTAKEN: // Strongly not taken lower bound
                    two_bit_ctr[pht_idx][entry_idx] = STRONG_NOTTAKEN;
                    break;
                case MAJOR_NOTTAKEN:
                case MINOR_NOTTAKEN:
                case WEAK_NOTTAKEN: 
                case WEAK_TAKEN:   
                case MINOR_TAKEN:   
                case MAJOR_TAKEN:   
                case STRONG_TAKEN: // Strongly taken upper bound
                    two_bit_ctr[pht_idx][entry_idx]--;
                    break;
                default:
                    // The first time update taken
                    two_bit_ctr[pht_idx][entry_idx] = WEAK_NOTTAKEN;
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
        unsigned int BH;
        unsigned int BA;

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
        unsigned int Get_Hist () {
                BH = bhsr_k_idx;
                return BH;
            };

        unsigned int Get_Branch (){
                BA = bi.address & ((1 << 16) - 1);
                return BA;
        };

};

////////////////////////////////////
// Selector
////////////////////////////////////

#define Components 2 //Branch Predictor Components
#define BHB 6 // Branch History Bits
#define BAB 4 // Branch Address Bits
#define Bindex 1<<(BHB+BAB) // VMT Branch index
#define Cindex 1<<Components // VMT Components index
#define Ctrbits 4 //Counter bits 3 bits is better
#define CtrMax (1<<Ctrbits)-1
#define CtrMin 0

int VMT [Bindex][Cindex] = {1<<(Ctrbits-1)};

class Vector_Mapping_Table{

    public:

    Vector_Mapping_Table(){}

    ~Vector_Mapping_Table(){}

    bool selection(unsigned long long bindex, unsigned long long cindex){
        return VMT[bindex][cindex]>>(Ctrbits-1);
    }

    void update_VMT(unsigned long long bindex, unsigned long long cindex, bool taken){
            int pred = VMT[bindex][cindex];
            
            if (taken){
                if(pred == CtrMax){
                    VMT[bindex][cindex] = CtrMax;
                }
                else{
                    VMT[bindex][cindex]++;
                }
            }
            else{
                if(pred == CtrMin){
                    VMT[bindex][cindex] = CtrMin;
                }
                else{
                    VMT[bindex][cindex]--;
                }
            }
    }
};


class selector{
    private:
        class Vector_Mapping_Table vmt;
        TAGE *tage_bp;
        twolevel_predictor *twolevel_bp;
        unsigned long long bindex,BHR;
        unsigned int cindex,PC;

    public:

        branch_update u;

        selector()
        {
            tage_bp = new TAGE();
	        twolevel_bp = new twolevel_predictor();
        }

        ~selector()
        {
            delete tage_bp;
            delete twolevel_bp;
        }

        branch_update *predict(branch_info &b)
        {
            bool bx = tage_bp->predict(b)->direction_prediction();
            bool cx = twolevel_bp->predict(b)->direction_prediction();
            BHR = ((twolevel_bp->Get_Hist())&((1<<BHB)-1))<<BAB;
            PC = (twolevel_bp->Get_Branch())&((1<<BAB)-1);
            bindex = BHR + PC;
            cindex = (bx<<(Components-1))+cx;
            u.direction_prediction(vmt.selection(bindex,cindex));
            u.target_prediction (0);
            return &u;
        }

        void update (branch_update * u, bool taken, unsigned int target)
        {
            // Here we will update both predictor
            tage_bp->update(u, taken, target);
            twolevel_bp->update(u, taken, target);
            vmt.update_VMT(bindex,cindex,taken);
        }
};
