#ifndef FRICTION_EXTRACTION_HPP
#define FRICTION_EXTRACTION_HPP

#include <vector>
#include <cmath>
#include "input_potential.hpp"

class friction_extraction
{
public:
        friction_extraction(std::vector<double> traj, double timestep, int bench_num, double beginning, double end)
                : traj_(traj),
                 timestep_(timestep),
                 bench_num_(bench_num),
                 beginning_(beginning),
                 end_(end)
        {}
        void forward_ffpt_calc();
        void forward_friction();



        double crossing_time(double boundary, int jump_index);

        std::vector<double> boundary_array;
        void state_generator(double beginning_, double end, int bench_num);
        void forward_fpt(int beginning_index);
        double forward_partition(int q_F);
        // q_F is the index of end point for partiion function integral.

        //




        void backward_fpt(int beginning_index);
        std::vector<std::vector<double>> ffpt;
        std::vector<std::vector<int>> ffpt_index;
        std::vector<bool> jump_state;
        std::vector<double> mffpt;
        std::vector<double> forward_friction_kernel;

        std::vector<double> afpt;
        std::vector<double> start_fpt;
        // fpt will record the totall time of FFPT
        // start_fpt records beginning time (beginning state) to this specific state.
        /* 
        std::vector<int> jump_time(bench_num_);
        std::vector<int> jump_cross_time(bench_num_);
        std::vector<double> mffpt(bench_num_);
        std::vector<double> mafpt(bench_num_);
        */ 
        void towards(double beginning, double end);


        // Save the discretize size; step = (begin - end)/num

private:
        std::vector<double> traj_;
        double timestep_;
        int bench_num_;
        int towards_;
        double end_;
        double beginning_;
        double step;



};
void friction_extraction::towards(double beginning, double end)
// Direction for forwards or backwards
{
        if(end>beginning)
                towards_  = 1;
        else
                towards_  = -1;

}
void friction_extraction::state_generator(double beginning, double end, int bench_num)
{
        double length = end-beginning;
        step   = length/bench_num;
        std::cout<<"Bench num is "<<step<<std::endl;
        for(int i=0; i<bench_num; i++)
        {
                boundary_array.push_back(beginning+step*i);
        }
        boundary_array.push_back(end_);
        ffpt.resize(bench_num);
        ffpt_index.resize(bench_num);
        jump_state.resize(bench_num);
}


double friction_extraction::crossing_time(double boundary, int  jump_index)
/*  *        the traj_ before jump_index 
 *   \
 *---------- boundary
 *     \
 *      *    the traj_[jump_index]
 * \ means the path, * means the data in traj_.
 * To get the crossing time, first calculate the velocity, then the time of jump_index
 * minus the time between the boundary and jump_index;
 */
{
        double distance_bound = fabs(traj_[jump_index]-boundary);
        double distance = distance_bound + fabs(boundary-traj_[jump_index-1]-boundary);
        double velocity = distance/timestep_;
        double time_difference = distance_bound/velocity;
        double jump_time = jump_index*timestep_;
        return jump_time - time_difference;
}

void friction_extraction::forward_fpt(int beginning_index)
{
        double beginning_time = crossing_time(beginning_, beginning_index);
        std::cout<<"this is the size"<<start_fpt.size()<<std::endl;
        std::cout<<"This is bench_num "<<bench_num_<<std::endl;
        for( int i=0; i<bench_num_; i++)
               start_fpt.push_back(beginning_time);
        int mid=1;
        bool above_down = true;
        double trans_time;
        double ffp_time;
        bool skip_status;
       // When true, means particle above boundary, when false, means below
        for(auto & b : boundary_array)
                std::cout<<"Boundary is "<<b<<std::endl;
        
        for(int i=beginning_index; i<traj_.size(); i++)
        {
                if(skip_status && traj_[i]>beginning_)
                        continue;
                else if(traj_[i]<beginning_&&traj_[i-1]>=beginning_)
                {
                        if(skip_status)
                        {
                                skip_status = false;
                        }
                        mid = bench_num_;
                        trans_time = crossing_time(boundary_array[0], i);
                        for(int j=0; j<jump_state.size(); j++)
                        {
                                if(jump_state[j]==false)
                                {
                                        jump_state[j] = true;
                                        start_fpt[j]  = trans_time;
                                }
                        }

                        // Initialize all the jump time;
                }
                if(traj_[i]>=boundary_array[mid-1])
                {
                        if(traj_[i]>boundary_array[mid])
                        {
                                // Update up and down state
                                if(traj_[i]>end_ &&jump_state[bench_num_-1]==false)
                                {
                                        skip_status = true;
                                        continue;
                                }


                                trans_time = crossing_time(boundary_array[mid], i);
                                if(jump_state[mid-1])
                                {
                                        ffpt[mid-1].push_back(trans_time-start_fpt[mid-1]);
                                        ffpt_index[mid-1].push_back(i);
                                        jump_state[mid-1] = false;
                                        // false means particle already passed the boundary
                                        // true means still waiting for count time.
                                        // After particle reach the state boundary. It is initialized.
                                }

                                if(mid<bench_num_)
                                        mid += 1;
                        }

                }
                else
                {
                        if (mid > 1)
                                mid -= 1;
                }
                
        }

}




void friction_extraction::forward_ffpt_calc()
{
        // Initializing the system
        towards(beginning_, end_);

        std::cout<<"Number is "<<bench_num_<<std::endl;
        state_generator(beginning_, end_, bench_num_);
        int index = 0;
        std::cout<<"Initializing"<<std::endl;
        while((traj_[index]-beginning_)*(traj_[index+1]-beginning_)>0)
        {
                index += 1;
        }
        std::cout<<"Find initial index"<<std::endl;
        double beginning_time = crossing_time(beginning_, index);
        std::cout<<"get beginning time"<<beginning_time<<std::endl;

        std::cout<<"before loop:"<<start_fpt.size()<<std::endl;
        // Fisrt First passage time and all first passage time;
        std::cout<<"FFPT calculating"<<std::endl;
        forward_fpt(index);
        double sum;
        std::cout<<"MFPT"<<std::endl;
        for(int i=0; i<bench_num_; i++)
        {
                for (auto & t : ffpt[i])
                        sum += t/ffpt[i].size();
                mffpt.push_back(sum);
                sum = 0;
        }
}
/*
 * Backward calculation
 */

void friction_extraction::backward_fpt(int beginning_index)
{
        double beginning_time = crossing_time(end_, beginning_index);
        std::cout<<"this is the size"<<start_fpt.size()<<std::endl;
        std::cout<<"This is bench_num "<<bench_num_<<std::endl;
        for( int i=0; i<bench_num_; i++)
               start_fpt.push_back(beginning_time);
        int mid=1;
        bool above_up = true;
        double trans_time;
        double ffp_time;
        bool skip_status=false;
       // When true, means particle below the upper boundary, when false, means above
        for(auto & b : boundary_array)
                std::cout<<"Boundary is "<<b<<std::endl;
        
        for(int i=beginning_index; i<traj_.size(); i++)
        {
                if(skip_status && traj_[i]<end_)
                        continue;
                // After a transition, everything just skipped, before it goes
                // back
                else if(traj_[i]>end_&&traj_[i-1]<=end_)
                {
                        if(skip_status)
                        {
                                skip_status = false;
                        }
                        mid = bench_num_;
                        trans_time = crossing_time(boundary_array[0], i);
                        for(int j=0; j<jump_state.size(); j++)
                        {
                                if(jump_state[j]==false)
                                {
                                        jump_state[j] = true;
                                        start_fpt[j]  = trans_time;
                                }
                        }

                        // Initialize all the jump time;
                }
                if(traj_[i]<=boundary_array[mid])
                {
                        if(traj_[i]<boundary_array[mid-1])
                        {
                                // Update up and down state
                                if(traj_[i]<beginning_ &&jump_state[0]==false)
                                {
                                        skip_status = true;
                                        continue;
                                }


                                trans_time = crossing_time(boundary_array[mid-1], i);
                                if(jump_state[mid])
                                {
                                        ffpt[mid].push_back(trans_time-start_fpt[mid-1]);
                                        ffpt_index[mid].push_back(i);
                                        jump_state[mid] = false;
                                        // false means particle already passed the boundary
                                        // true means still waiting for count time.
                                        // After particle reach the state boundary. It is initialized.
                                }

                                if(0<mid)
                                        mid -= 1;
                        }

                }
                else
                {
                        if (mid<bench_num_)
                                mid += 1;
                }
                
        }

}




void friction_extraction::forward_ffpt_calc()
{
        // Initializing the system
        towards(beginning_, end_);

        std::cout<<"Number is "<<bench_num_<<std::endl;
        state_generator(beginning_, end_, bench_num_);
        int index = 0;
        std::cout<<"Initializing"<<std::endl;
        while((traj_[index]-beginning_)*(traj_[index+1]-beginning_)>0)
        {
                index += 1;
        }
        std::cout<<"Find initial index"<<std::endl;
        double beginning_time = crossing_time(beginning_, index);
        std::cout<<"get beginning time"<<beginning_time<<std::endl;

        std::cout<<"before loop:"<<start_fpt.size()<<std::endl;
        // Fisrt First passage time and all first passage time;
        std::cout<<"FFPT calculating"<<std::endl;
        forward_fpt(index);
        double sum;
        std::cout<<"MFPT"<<std::endl;
        for(int i=0; i<bench_num_; i++)
        {
                for (auto & t : ffpt[i])
                        sum += t/ffpt[i].size();
                mffpt.push_back(sum);
                sum = 0;
        }
}



// Calculate the partition integral Z1 =∫^qF_qmin dq e^(βU(q))
// rectangle rule

double Boltzmann(double x)
{
        double beta = 1;
        //beta = 1/k_BT
        return exp(-beta*input_potential(x));
}
void friction_extraction::forward_friction()
{
// Calculate the partition integral Z1 =∫^qF_qmin dq e^(βU(q))
// k_BT = 1
        std::vector<double> Z1;
        double integral;
        double friction;
        
        for(int i=0; i<bench_num_/2; i++)
                integral += 0.5*step*(Boltzmann(boundary_array[0]-step*i)+Boltzmann(boundary_array[1]-step*i));
        friction = Boltzmann(boundary_array[1])/integral;
        friction = friction*mffpt[0]/step;
        forward_friction_kernel.push_back(friction); 
        // partial derivation.
 
        for(int i=1; i<bench_num_; i++)
        {
                integral += 0.5*step*(Boltzmann(boundary_array[i])+Boltzmann(boundary_array[i+1]));
                friction = Boltzmann(boundary_array[i+1])/integral;
                friction = friction*(mffpt[i]-mffpt[i-1])/step;
                forward_friction_kernel.push_back(friction);
                Z1.push_back(integral);
        }
        std::cout<<"This is size of friction: "<<forward_friction_kernel.size()<<std::endl;
}






#endif
