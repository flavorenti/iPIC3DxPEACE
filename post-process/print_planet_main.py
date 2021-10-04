#### rappel of the routines ####
from print_planet_routines import *

# main has the arguments: 1) Production run number (0,1,2,0-bigbox)
#                         2) factor proportional to filter width and scaling of the box (PR0,1,2 -> fact=1,2,4)
#                         2) nrun
#                         3) Ncycle start
#                         4) Ncycle end
#                         5) Ncycle step 
#                         6) Interval pictures (default 1-6)

args=[]

# PR0 run0
start = 0
end   = 5000#3700
for i in range(start,end,100): 
    args.append(['0',1,'_new',i,i+100,100,1,7])

'''
# PR0 run1
start = 3800
end   = 5700
for i in range(start,end,100):
    args.append(['0',1,'1',i,i+100,100,1,6])
'''
'''
# PR1 run0
start = 0
end   = 5100
for i in range(start,end,100):
    args.append(['1',2,0,i,i+100,100,1,6])
# PR1 run1
start = 5100
end   = 8500
for i in range(start,end,100):
    args.append(['1',2,1,i,i+100,100,1,6])
'''
'''
# PR2 run0 fig1
#args.append(['2',4,0,0,1000,100,1,1])
# PR2 run0 fig2-3
args.append(['2',4,0,0,1100,100,2,2])
# PR2 run0 fig4-5
#args.append(['2',4,0,0,1000,100,4,5])
# PR2 run1 fig2-3
args.append(['2',4,1,1100,1500,100,2,2])
# PR2 run2-13 
start = 1500
end = 2000
for irun in range(2,14):
    #args.append(['2',4,irun,start,end,100,1,1])
    args.append(['2',4,irun,start,end,100,2,2])
    #args.append(['2',4,irun,start,end,100,4,5])
    start+=500
    end+=500
# PR2 run14-27
start = 7500
end = 7750
for irun in range(14,28):
    #args.append(['2',4,irun,start+50*(irun%2),end,100,1,1]) 
    args.append(['2',4,irun,start+50*(irun%2),end,100,2,3])
    #args.append(['2',4,irun,start+50*(irun%2),end,100,4,5])
    start+=250
    end+=250
'''
inow=40
main(args[inow])

# for PR0-bigbox
#main('0-bigbox',2,0,100,5100,100)
#main('0-bigbox',2,1,5100,10000,100)
