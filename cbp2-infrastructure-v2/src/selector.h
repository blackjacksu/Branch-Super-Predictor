#include "my_predictor.h"
#include "twolevel_predictor.h"

#define Components 2 //Branch Predictor Components
#define BHB 6
#define BAB 4
#define Bindex 1<<(BHB+BAB) // VMT Branch index
#define Cindex 1<<Components // VMT Components index
#define Ctrbits 2 //Counter bits
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
