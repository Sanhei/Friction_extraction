# Friction_extraction
Depend on the Fokker Plank equation, we can extract the position dependent friction from the first first passage time
Only for internal using. (AG Netz)

1. Changing the input file in main.cpp. Find the line
```
inputfile.open("0.00100.txt")
```
change 0.00100.txt to the file name you use.

2. find the line 
```
int bench_num=50;
```
This is how many bins you want to discretize, e.g. from -1 to 1 and ben_num=50, it will calculate the First-First passage time (FFPT) from
-1 to $-1+\frac{1-(-1)}{50}*i$, i=1, 2, 3. You may need to try this several times to get proper FFPT for the later friction calculation.

3. find the line
```
friction_extraction FE(traj, 0.0001, bench_num, -1, 1);
```
Change 0.0001 to your time lag.
Change -1 to your start position.
Change 1 to your end position.

4. Compile.
After saving the files, in the terminal do
```
bash run.sh
```
Then, if you are not using double well potential, then you only need to check the ffpt.txt, which record your FFPT. And in the python you may create a corresponding position array np.linspace(start, end, (end-start)/ben_num).





