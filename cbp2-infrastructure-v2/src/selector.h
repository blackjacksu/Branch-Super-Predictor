#include "my_predictor.h"
#include "twolevel_predictor.h"


class selector{
    private:
        branch_predictor *tage_bp;
        branch_predictor *twolevel_bp;

    
    public: 
        selector()
        {
            tage_bp = new my_predictor();
	        twolevel_bp = new twolevel_predictor();
        }

        ~selector()
        {
            delete tage_bp;
            delete twolevel_bp;
        }

        branch_update *predict(branch_info &b)
        {
            branch_update *bx;
            if (0) // Selector logic
            {
                // Select tage predicotr
                bx = tage_bp->predict(b);
            }
            else
            {
                // Select two level predictor result
                bx = twolevel_bp->predict(b);
            }
            return bx;
        }

        void update (branch_update * u, bool taken, unsigned int target)
        {
            // Here we will update both predictor
            tage_bp->update(u, taken, target);
            twolevel_bp->update(u, taken, target);
        }

        // 
};