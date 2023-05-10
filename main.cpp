#include <fstream>
#include <iostream>

#include <vector>
#include "friction_extraction.hpp"
int main()
{
        std::vector<double> traj;
        std::ifstream inputfile;
        inputfile.open("0.001000.txt");
        double x;
        int bench_num = 50;
        while(inputfile>>x)
                traj.push_back(x);
        friction_extraction FE(traj, 0.0001, bench_num, -1, 1);
        FE.forward_ffpt_calc();
        FE.forward_friction();

        
        std::ofstream ffpt;
        ffpt.open("ffpt.txt");
        for(int i=bench_num-1; i<FE.ffpt.size(); i++)
        {
                for(auto & t : FE.ffpt[i])
                        ffpt<<t<<" ";
                ffpt<<std::endl;
        }
        ffpt.close();

        std::ofstream mffpt;
        mffpt.open("mffpt.txt");
        for(int i=0; i<FE.mffpt.size(); i++)
        {
                mffpt<<FE.mffpt[i]<<std::endl;
        }
        mffpt.close();

        std::ofstream ffk;
        ffk.open("forward_friction_kernel.txt");
        for(int i=0; i<FE.forward_friction_kernel.size(); i++)
        {
                ffk<<FE.forward_friction_kernel[i]<<std::endl;
        }
        ffk.close();



}
