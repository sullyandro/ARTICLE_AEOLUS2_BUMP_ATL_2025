"""Defines linear equations of 2.0 layer of Aeolus2.0, rostami@pik-potsdam.de"""
# This script sets up linear equation of mcTRSW dynamical core

import numpy             as np
import scipy.sparse      as sparse
#import matplotlib.pyplot as plt
#from scipy.linalg         import eig
#from mpl_toolkits.basemap import Basemap

# parameters:
# g: gravity
# H: equilibrium height
# Om: rotation rate

# equations in 'dry' case (km= grad_m  h1 is practically eta1, as is written in the model's manuscript):
# dt(u1m) + B1*km*h1 + B1*km*h2 + B1*km*h3 +0.5*H1*km*b'1 - 2*Om*i*C*u1m = - (u1.grad u1)_m +0.5*h1. km*b'1 - km*(h1*b'1) - b'1. km* h3 - b'1.km*h2     !
# dt(u1p) + B1*kp*h1 + B1*kp*h2 + B1*kp*h3 +0.5*H1*kp*b'1 + 2*Om*i*C*u1p = - (u1.grad u1)_p +0.5*h1. kp*b'1 - kp*(h1*b'1) - b'1. kp*h3  - b'1.kp*h2
# dt(h1)  + H*kp*u1m + H*km*u1p                                          = - kp*(h*u1)_m - km*(h*u1)_p

# dt(u2m) + B1*km*h1 + B2*km*h2 + B2*km*h3 + H1*kp*b'1 + 0.5*H2*kp*b'2 - 2*Om*i*C*u2m = - (u2.grad u2)_m + 0.5*h2*km*b'2 - km*(h2*b'2) - km*(h1*b'1) - b'2*km*h3
# dt(u2p) + B1*km*h1 + B2*km*h2 + B2*km*h3 + H1*kp*b'1 + 0.5*H2*kp*b'2 + 2*Om*i*C*u2p = - (u2.grad u2)_p + 0.5*h2*kp*b'2 - kp*(h2*b'2) - kp*(h1*b'1) - b'2*kp*h3
# dt(h) + H*kp*um + H*km*up = - kp*(h*u)_m - km*(h*u)_p
# dt(b1) = - (u1.grad b1)
#-------------------------------------------------------------------------------

# variable orders: um1 up1 h1 um2 up2 h2 q1 b1 b2 q2 w1 w2 
#                  0   1   2  3   4   5  6  7  8  9  10 11

def shallow_water(S,m,params):
    """Defines M, L matrices for advection"""
        
#    g,H,Om,a  = params[0],params[1],params[2],params[3]
#    Mm,Lm = eq.shallow_water(S,m,[g0,B2,B1,H1,H2,Om,a,nu,kappa,Q01])
    g0,B3,B2,B1,H1,H2,Om,a,nu,kappa,Q01,Q02  = params[0],params[1],params[2],params[3],params[4],params[5],params[6],params[7],params[8],params[9],params[10],params[11]
    g=g0*np.array([[B1, B1, B1],[B1, B2, B2],[B1,B2,B3]])

#  layer 1 
    # (u1-,u1-)
    L00   = -2*Om*1j*S.op('C',m,-1)-nu*((S.op('k-',m,0)).dot(S.op('k+',m,-1))+(S.op('k+',m,-2)).dot(S.op('k-',m,-1)))
    # (u1-,u1+)
    L01   = S.zeros(m,-1,1)
    # (u1-,h1)
    L02   = g[0][0]*S.op('k-',m,0)/a
    # (u1-,u2-)
    L03   = S.zeros(m,-1,-1)
    # (u1-,u2+)
    L04   = S.zeros(m,-1,1)
    # (u1-,h2)
    L05   = g[0][1]*S.op('k-',m,0)/a
    # (u1-,q1)
    L06   = S.zeros(m,-1,0)
    # (u1-,b1)
    L07   = 0.5*H1*S.op('k-',m,0)/a # S.zeros(m,-1,0) # 
    # (u1-,b2)    
    L08   = S.zeros(m,-1,0) 
    # (u1-,q2)
    L09   = S.zeros(m,-1,0)    
    # (u1-,w1)
    L0_10 = S.zeros(m,-1,0)    
    # (u1-,w2)
    L0_11 = S.zeros(m,-1,0) 
           
    # (u1+,u1-)
    L10   = S.zeros(m,1,-1)
    # (u1+,u1+)
    L11   = 2*Om*1j*S.op('C',m,+1)-nu*((S.op('k-',m,2)).dot(S.op('k+',m,1))+(S.op('k+',m,0)).dot(S.op('k-',m,1)))
    # (u1+,h1)
    L12   = g[0][0]*S.op('k+',m,0)/a
    # (u1+,u2-)
    L13   = S.zeros(m,1,-1)
    # (u1+,u2+)
    L14   = S.zeros(m,1,+1)
    # (u1+,h2)
    L15   = g[0][1]*S.op('k+',m,0)/a
    # (u1+,q1)
    L16   = S.zeros(m,1,0)
    # (u1+,b1)
    L17   = 0.5*H1*S.op('k+',m,0)/a # S.zeros(m,1,0) # 
    # (u1+,b2)    
    L18   = S.zeros(m,-1,0)
    # (u1+,q2)
    L19   = S.zeros(m,1,0)    
    # (u1+,w1)
    L1_10 = S.zeros(m,1,0)
    # (u1+,w2)
    L1_11 = S.zeros(m,1,0)  
    
    # (h1,u1-)
    L20   = H1*S.op('k+',m,-1)/a
    # (h1,u1+)
    L21   = H1*S.op('k-',m,1)/a
    # (h1,h1)
    L22   = -kappa*((S.op('k-',m,1)).dot(S.op('k+',m,0))+(S.op('k+',m,-1)).dot(S.op('k-',m,0)))
    # (h1,u2-)
    L23   = S.zeros(m,0,-1)
    # (h1,u2+)
    L24   = S.zeros(m,0,+1)
    # (h1,h2)
    L25   = S.zeros(m,0,0)
    # (h1,q1)
    L26   = S.zeros(m,0,0)
   # (h1,b1)
    L27   = S.zeros(m,0,0)
   # (h1,b2)    
    L28   = S.zeros(m,0,0) 
    # (h1,q2)
    L29   = S.zeros(m,0,0)
    # (h1,w1)
    L2_10 = S.zeros(m,0,0)
    # (h1,w2)
    L2_11 = S.zeros(m,0,0)
        
#  layer 2 -------------------------------------------
    # (u2-,u1-)
    L30   = S.zeros(m,-1,-1)
    # (u2-,u1+)
    L31   = S.zeros(m,-1,1)
    # (u2-,h1)
    L32   = g[1][0]*S.op('k-',m,0)/a
    # (u2-,u2-)
    L33   = -2*Om*1j*S.op('C',m,-1)-nu*((S.op('k-',m,0)).dot(S.op('k+',m,-1))+(S.op('k+',m,-2)).dot(S.op('k-',m,-1)))
    # (u2-,u2+)
    L34   = S.zeros(m,-1,1)
    # (u2-,h2)
    L35   = g[1][1]*S.op('k-',m,0)/a
    # (u2-,q1)
    L36   = S.zeros(m,-1,0)
    # (u2-,b1)
    L37   = H1*S.op('k-',m,0)/a 
    # (u2-,b2)
    L38   = 0.5*H2*S.op('k-',m,0)/a 
    # (u2-,q2)
    L39   = S.zeros(m,-1,0)
    # (u2-,w1)
    L3_10   = S.zeros(m,-1,0)
    # (u2-,w2)
    L3_11 = S.zeros(m,-1,0)

    # (u2+,u1-)
    L40   = S.zeros(m,1,-1)
    # (u2+,u1+)
    L41   = S.zeros(m,1,1)
    # (u2+,h1)
    L42   = g[1][0]*S.op('k+',m,0)/a
    # (u2+,u2-)
    L43   = S.zeros(m,+1,-1)
    # (u2+,u2+)
    L44   = 2*Om*1j*S.op('C',m,+1)-nu*((S.op('k-',m,2)).dot(S.op('k+',m,1))+(S.op('k+',m,0)).dot(S.op('k-',m,1)))
    # (u2+,h2)
    L45   = g[1][1]*S.op('k+',m,0)/a
    # (u2+,q1)
    L46   = S.zeros(m,1,0)
    # (u2+,b1)
    L47   = H1*S.op('k+',m,0)/a 
    # (u2+,b2)
    L48   = 0.5*H2*S.op('k+',m,0)/a
    # (u2+,q2)
    L49   = S.zeros(m,1,0)
    # (u2+,w1)
    L4_10 = S.zeros(m,1,0)
    # (u2+,w2)
    L4_11 = S.zeros(m,1,0)
    
    # (h2,u1-)
    L50   = S.zeros(m,0,-1)
    # (h2,u1+)
    L51   = S.zeros(m,0,+1)
    # (h2,h1)
    L52   = S.zeros(m,0,0)
    # (h2,u2-)
    L53   = H2*S.op('k+',m,-1)/a
    # (h2,u2+)
    L54   = H2*S.op('k-',m,1)/a
    # (h2,h2)
    L55   = -kappa*((S.op('k-',m,1)).dot(S.op('k+',m,0))+(S.op('k+',m,-1)).dot(S.op('k-',m,0))) #S.zeros(m,0,0)
    # (h2,q1)
    L56   = S.zeros(m,0,0)
    # (h2,b1)
    L57   = S.zeros(m,0,0)
    # (h2,b2)
    L58   = S.zeros(m,0,0)
    # (h2,q2)
    L59   = S.zeros(m,0,0)
    # (h2,w1)
    L5_10 = S.zeros(m,0,0)
    # (h2,w2)
    L5_11 = S.zeros(m,0,0)
        
    # (q1,u1-)
    L60   = Q01*S.op('k+',m,-1)/a
    # (q1,u1+)
    L61   = Q01*S.op('k-',m,1)/a
    # (q1,h1)
    L62   = S.zeros(m,0,0)
    # (q1,u2-)
    L63   = S.zeros(m,0,-1)
    # (q1,u2+)
    L64   = S.zeros(m,0,+1)
    # (q1,h2)
    L65   = S.zeros(m,0,0)
    # (q1,q1)
    L66   = -nu*((S.op('k-',m,1)).dot(S.op('k+',m,0))+(S.op('k+',m,-1)).dot(S.op('k-',m,0)))
    # (q1,b1)
    L67   = S.zeros(m,0,0)    
    # (q1,b2)
    L68   = S.zeros(m,0,0)
    # (q1,q2)
    L69   = S.zeros(m,0,0)    
    # (q1,w1)
    L6_10 = S.zeros(m,0,0)
    # (q1,w2)
    L6_11 = S.zeros(m,0,0)

    # (b1,u1-)
    L70   = S.zeros(m,0,-1)
    # (b1,u1+)
    L71   = S.zeros(m,0,+1)
    # (b1,h1)
    L72   = S.zeros(m,0,0)
    # (b1,u2-)
    L73   = S.zeros(m,0,-1)
    # (b1,u2+)
    L74   = S.zeros(m,0,+1)
    # (b1,h2)
    L75   = S.zeros(m,0,0)
    # (b1,q1)
    L76   = S.zeros(m,0,0)
    # (b1,b1)
    L77   = -nu*((S.op('k-',m,1)).dot(S.op('k+',m,0))+(S.op('k+',m,-1)).dot(S.op('k-',m,0)))    #S.zeros(m,0,0) 
    # (b1,b2) 
    L78   = S.zeros(m,0,0) 
    # (b1,q2) 
    L79   = S.zeros(m,0,0)     
    # (b1,w1)
    L7_10 = S.zeros(m,0,0) 
    # (b1,w2) 
    L7_11 = S.zeros(m,0,0) 
    
    # (b2,u1-)
    L80   = S.zeros(m,0,-1)
    # (b2,u1+)
    L81   = S.zeros(m,0,+1)
    # (b2,h1)
    L82   = S.zeros(m,0,0)
    # (b2,u2-)
    L83   = S.zeros(m,0,-1)
    # (b2,u2+)
    L84   = S.zeros(m,0,+1)
    # (b2,h2)
    L85   = S.zeros(m,0,0)
    # (b2,q1)
    L86   = S.zeros(m,0,0)
    # (b2,b1)
    L87   = S.zeros(m,0,0) 
    # (b2,b2) 
    L88   = -nu*((S.op('k-',m,1)).dot(S.op('k+',m,0))+(S.op('k+',m,-1)).dot(S.op('k-',m,0)))    #S.zeros(m,0,0) 
    # (b2,q2) 
    L89   = S.zeros(m,0,0) 
    # (b2,w1)
    L8_10 = S.zeros(m,0,0) 
    # (b2,w2) 
    L8_11 = S.zeros(m,0,0) 

    # (q2,u1-)     
    L90   = Q01*S.op('k+',m,-1)/a
    # (q2,u1+)
    L91   = Q01*S.op('k-',m,1)/a
    # (q2,h1)
    L92   = S.zeros(m,0,0)
    # (q2,u2-)
    L93   = S.zeros(m,0,-1)
    # (q2,u2+)
    L94   = S.zeros(m,0,+1)
    # (q2,h2)
    L95   = S.zeros(m,0,0)
    # (q2,q1)
    L96   = S.zeros(m,0,0)
    # (q2,b1)
    L97   = S.zeros(m,0,0) 
    # (q2,b2) 
    L98   = S.zeros(m,0,0) 
    # (q2,q2) 
    L99   = -nu*((S.op('k-',m,1)).dot(S.op('k+',m,0))+(S.op('k+',m,-1)).dot(S.op('k-',m,0)))
    # (q2,w1)
    L9_10 = S.zeros(m,0,0)
    # (q2,w2)
    L9_11 = S.zeros(m,0,0)   
    
    # (w1,u1-)
    L10_0   = S.zeros(m,0,-1)
    # (w1,u1+)
    L10_1   = S.zeros(m,0,+1)
    # (w1,h1)
    L10_2   = S.zeros(m,0,0)
    # (w1,u2-)
    L10_3   = S.zeros(m,0,-1)
    # (w1,u2+)
    L10_4   = S.zeros(m,0,+1)
    # (w1,h2)
    L10_5   = S.zeros(m,0,0)
    # (w1,q1)
    L10_6   = S.zeros(m,0,0)
    # (w1,b1)
    L10_7   = S.zeros(m,0,0) 
    # (w1,b2) 
    L10_8   = S.zeros(m,0,0) 
    # (w1,q2)
    L10_9   = S.zeros(m,0,0) 
    # (w1,w1) 
    L10_10  = -nu*((S.op('k-',m,1)).dot(S.op('k+',m,0))+(S.op('k+',m,-1)).dot(S.op('k-',m,0)))#S.zeros(m,0,0)       
    # (w1,w2) 
    L10_11  = S.zeros(m,0,0)

    # (w2,u1-)
    L11_0   = S.zeros(m,0,-1)
    # (w2,u1+)
    L11_1   = S.zeros(m,0,+1)
    # (w2,h1)
    L11_2   = S.zeros(m,0,0)
    # (w2,u2-)
    L11_3   = S.zeros(m,0,-1)
    # (w2,u2+)
    L11_4   = S.zeros(m,0,+1)
    # (w2,h2)
    L11_5   = S.zeros(m,0,0)
    # (w2,q1)
    L11_6   = S.zeros(m,0,0)
    # (w2,b1)
    L11_7   = S.zeros(m,0,0) 
    # (w2,b2) 
    L11_8   = S.zeros(m,0,0) 
    # (w2,q2)
    L11_9   = S.zeros(m,0,0) 
    # (w2,w1) 
    L11_10  = S.zeros(m,0,0)       
    # (w2,w2) 
    L11_11  = -nu*((S.op('k-',m,1)).dot(S.op('k+',m,0))+(S.op('k+',m,-1)).dot(S.op('k-',m,0)))#S.zeros(m,0,0)
    
    L = sparse.bmat([[L00,  L01,   L02,   L03,   L04,   L05,   L06,   L07,   L08,   L09,   L0_10,  L0_11], 
                    [L10,   L11,   L12,   L13,   L14,   L15,   L16,   L17,   L18,   L19,   L1_10,  L1_11],              
                    [L20,   L21,   L22,   L23,   L24,   L25,   L26,   L27,   L28,   L29,   L2_10,  L2_11], 
                    [L30,   L31,   L32,   L33,   L34,   L35,   L36,   L37,   L38,   L39,   L3_10,  L3_11], 
                    [L40,   L41,   L42,   L43,   L44,   L45,   L46,   L47,   L48,   L49,   L4_10,  L4_11], 
                    [L50,   L51,   L52,   L53,   L54,   L55,   L56,   L57,   L58,   L59,   L5_10,  L5_11],  
                    [L60,   L61,   L62,   L63,   L64,   L65,   L66,   L67,   L68,   L69,   L6_10,  L6_11], 
                    [L70,   L71,   L72,   L73,   L74,   L75,   L76,   L77,   L78,   L79,   L7_10,  L7_11], 
                    [L80,   L81,   L82,   L83,   L84,   L85,   L86,   L87,   L88,   L89,   L8_10,  L8_11], 
                    [L90,   L91,   L92,   L93,   L94,   L95,   L96,   L97,   L98,   L99,   L9_10,  L9_11],
                    [L10_0, L10_1, L10_2, L10_3, L10_4, L10_5, L10_6, L10_7, L10_8, L10_9, L10_10, L10_11],
                    [L11_0, L11_1, L11_2, L11_3, L11_4, L11_5, L11_6, L11_7, L11_8, L11_9, L11_10, L11_11]])

    # u-
    R00   = S.op('I',m,-1)
    R01   = S.zeros(m,-1,1)
    R02   = S.zeros(m,-1,0)
    R03   = S.zeros(m,-1,-1)
    R04   = S.zeros(m,-1,1)
    R05   = S.zeros(m,-1,0)
    R06   = S.zeros(m,-1,0)
    R07   = S.zeros(m,-1,0)
    R08   = S.zeros(m,-1,0)
    R09   = S.zeros(m,-1,0)   
    R0_10 = S.zeros(m,-1,0)   
    R0_11 = S.zeros(m,-1,0)   
             
    # u+
    R10   = S.zeros(m,1,-1)
    R11   = S.op('I',m,1)
    R12   = S.zeros(m,1,0)
    R13   = S.zeros(m,1,-1)
    R14   = S.zeros(m,1,1)
    R15   = S.zeros(m,1,0)
    R16   = S.zeros(m,1,0)
    R17   = S.zeros(m,1,0)
    R18   = S.zeros(m,1,0)
    R19   = S.zeros(m,1,0)  
    R1_10 = S.zeros(m,1,0)  
    R1_11 = S.zeros(m,1,0)  
             
    # h1
    R20   = S.zeros(m,0,-1)
    R21   = S.zeros(m,0,1)
    R22   = S.op('I',m,0)
    R23   = S.zeros(m,0,-1)
    R24   = S.zeros(m,0,1)
    R25   = S.zeros(m,0,0)
    R26   = S.zeros(m,0,0)
    R27   = S.zeros(m,0,0)
    R28   = S.zeros(m,0,0)
    R29   = S.zeros(m,0,0)
    R2_10 = S.zeros(m,0,0)
    R2_11 = S.zeros(m,0,0)    
            
    # u-
    R30   = S.zeros(m,-1,-1)
    R31   = S.zeros(m,-1,1)
    R32   = S.zeros(m,-1,0)
    R33   = S.op('I',m,-1)
    R34   = S.zeros(m,-1,1)
    R35   = S.zeros(m,-1,0)
    R36   = S.zeros(m,-1,0)
    R37   = S.zeros(m,-1,0)
    R38   = S.zeros(m,-1,0)
    R39   = S.zeros(m,-1,0)
    R3_10 = S.zeros(m,-1,0)
    R3_11 = S.zeros(m,-1,0)
                
    # u+
    R40   = S.zeros(m,1,-1)
    R41   = S.zeros(m,1,1)
    R42   = S.zeros(m,1,0)
    R43   = S.zeros(m,1,-1)
    R44   = S.op('I',m,1)
    R45   = S.zeros(m,1,0)
    R46   = S.zeros(m,1,0)
    R47   = S.zeros(m,1,0)
    R48   = S.zeros(m,1,0)
    R49   = S.zeros(m,1,0)
    R4_10 = S.zeros(m,1,0)
    R4_11 = S.zeros(m,1,0)    
            
    # h2
    R50   = S.zeros(m,0,-1)
    R51   = S.zeros(m,0,1)
    R52   = S.zeros(m,0,0)
    R53   = S.zeros(m,0,-1)
    R54   = S.zeros(m,0,1)
    R55   = S.op('I',m,0)
    R56   = S.zeros(m,0,0)
    R57   = S.zeros(m,0,0)
    R58   = S.zeros(m,0,0)
    R59   = S.zeros(m,0,0)  
    R5_10 = S.zeros(m,0,0)    
    R5_11 = S.zeros(m,0,0)        
        
    # q1
    R60   = S.zeros(m,0,-1)
    R61   = S.zeros(m,0,1)
    R62   = S.zeros(m,0,0)
    R63   = S.zeros(m,0,-1)
    R64   = S.zeros(m,0,1)
    R65   = S.zeros(m,0,0)
    R66   = S.op('I',m,0)
    R67   = S.zeros(m,0,0)
    R68   = S.zeros(m,0,0)
    R69   = S.zeros(m,0,0)   
    R6_10 = S.zeros(m,0,0)    
    R6_11 = S.zeros(m,0,0)        
        
    # b1
    R70   = S.zeros(m,0,-1)
    R71   = S.zeros(m,0,1)
    R72   = S.zeros(m,0,0)
    R73   = S.zeros(m,0,-1)
    R74   = S.zeros(m,0,1)
    R75   = S.zeros(m,0,0)
    R76   = S.zeros(m,0,0)
    R77   = S.op('I',m,0)
    R78   = S.zeros(m,0,0)
    R79   = S.zeros(m,0,0)
    R7_10 = S.zeros(m,0,0)
    R7_11 = S.zeros(m,0,0)    
                
    # b2
    R80   = S.zeros(m,0,-1)
    R81   = S.zeros(m,0,1)
    R82   = S.zeros(m,0,0)
    R83   = S.zeros(m,0,-1)
    R84   = S.zeros(m,0,1)
    R85   = S.zeros(m,0,0)
    R86   = S.zeros(m,0,0)
    R87   = S.zeros(m,0,0)
    R88   = S.op('I',m,0)
    R89   = S.zeros(m,0,0)
    R8_10 = S.zeros(m,0,0)
    R8_11 = S.zeros(m,0,0)
    
    # q2
    R90   = S.zeros(m,0,-1)
    R91   = S.zeros(m,0,1)
    R92   = S.zeros(m,0,0)
    R93   = S.zeros(m,0,-1)
    R94   = S.zeros(m,0,1)
    R95   = S.zeros(m,0,0)
    R96   = S.zeros(m,0,0)
    R97   = S.zeros(m,0,0)
    R98   = S.zeros(m,0,0)
    R99   = S.op('I',m,0) 
    R9_10 = S.zeros(m,0,0)   
    R9_11 = S.zeros(m,0,0)     
    
    # w1
    R10_0   = S.zeros(m,0,-1)
    R10_1   = S.zeros(m,0,1)
    R10_2   = S.zeros(m,0,0)
    R10_3   = S.zeros(m,0,-1)
    R10_4   = S.zeros(m,0,1)
    R10_5   = S.zeros(m,0,0)
    R10_6   = S.zeros(m,0,0)
    R10_7   = S.zeros(m,0,0)
    R10_8   = S.zeros(m,0,0)
    R10_9   = S.zeros(m,0,0)                   
    R10_10  = S.op('I',m,0) 
    R10_11  = S.zeros(m,0,0)                   

    # w2
    R11_0   = S.zeros(m,0,-1)
    R11_1   = S.zeros(m,0,1)
    R11_2   = S.zeros(m,0,0)
    R11_3   = S.zeros(m,0,-1)
    R11_4   = S.zeros(m,0,1)
    R11_5   = S.zeros(m,0,0)
    R11_6   = S.zeros(m,0,0)
    R11_7   = S.zeros(m,0,0)
    R11_8   = S.zeros(m,0,0)
    R11_9   = S.zeros(m,0,0)                   
    R11_10  = S.zeros(m,0,0)   
    R11_11  = S.op('I',m,0) 
                
    R = sparse.bmat([  [R00,    R01,    R02,    R03,   R04,    R05,    R06,    R07,    R08,    R09,    R0_10,    R0_11],
                       [R10,    R11,    R12,    R13,   R14,    R15,    R16,    R17,    R18,    R19,    R1_10,    R1_11], 
                       [R20,    R21,    R22,    R23,   R24,    R25,    R26,    R27,    R28,    R29,    R2_10,    R2_11], 
                       [R30,    R31,    R32,    R33,   R34,    R35,    R36,    R37,    R38,    R39,    R3_10,    R3_11], 
                       [R40,    R41,    R42,    R43,   R44,    R45,    R46,    R47,    R48,    R49,    R4_10,    R4_11], 
                       [R50,    R51,    R52,    R53,   R54,    R55,    R56,    R57,    R58,    R59,    R5_10,    R5_11], 
                       [R60,    R61,    R62,    R63,   R64,    R65,    R66,    R67,    R68,    R69,    R6_10,    R6_11], 
                       [R70,    R71,    R72,    R73,   R74,    R75,    R76,    R77,    R78,    R79,    R7_10,    R7_11], 
                       [R80,    R81,    R82,    R83,   R84,    R85,    R86,    R87,    R88,    R89,    R8_10,    R8_11], 
                       [R90,    R91,    R92,    R93,   R94,    R95,    R96,    R97,    R98,    R99,    R9_10,    R9_11],
                       [R10_0,  R10_1,  R10_2,  R10_3, R10_4,  R10_5,  R10_6,  R10_7,  R10_8,  R10_9,  R10_10,   R10_11],
                       [R11_0,  R11_1,  R11_2,  R11_3, R11_4,  R11_5,  R11_6,  R11_7,  R11_8,  R11_9,  R11_10,   R11_11]])

    return R,L

class StateVector:

    def __init__(self,u1,h1,u2,h2,q1,b1,b2,q2,w1,w2):
        self.data  = []
        
        self.m_min = u1.S.m_min
        self.m_max = u1.S.m_max
        
        self.len_u1 = []
        self.len_h1 = []
        self.len_u2 = []
        self.len_h2 = []
        self.len_q1 = []
        self.len_b1 = []       
        self.len_b2 = [] 
        self.len_q2 = [] 
        self.len_w1 = []  # *** up to 'one' to the last variable 
            
        for m in range(self.m_min,self.m_max+1):
            m_local = m - self.m_min
            self.len_u1.append(u1['c'][m_local].shape[0])
            self.len_h1.append(h1['c'][m_local].shape[0])
            self.len_u2.append(u2['c'][m_local].shape[0])
            self.len_h2.append(h2['c'][m_local].shape[0])
            self.len_q1.append(q1['c'][m_local].shape[0])
            self.len_b1.append(b1['c'][m_local].shape[0])  
            self.len_b2.append(b2['c'][m_local].shape[0]) 
            self.len_q2.append(q2['c'][m_local].shape[0])
            self.len_w1.append(w1['c'][m_local].shape[0])  # *** up to 'one' to the last variable 
                                              
            self.data.append(np.concatenate((u1['c'][m_local], h1['c'][m_local], u2['c'][m_local], h2['c'][m_local], q1['c'][m_local], b1['c'][m_local], b2['c'][m_local], q2['c'][m_local], w1['c'][m_local], w2['c'][m_local])))

    def pack(self,u1,h1,u2,h2,q1,b1,b2,q2,w1,w2):
        for m in range(self.m_min,self.m_max+1):
            m_local            = m - self.m_min
            self.data[m_local] = np.concatenate((u1['c'][m_local], h1['c'][m_local], u2['c'][m_local], h2['c'][m_local], q1['c'][m_local], b1['c'][m_local],  b2['c'][m_local], q2['c'][m_local],  w1['c'][m_local],  w2['c'][m_local]))

    def unpack(self,u1,h1,u2,h2,q1,b1,b2,q2,w1,w2):
        u1.layout='c'
        h1.layout='c'
        u2.layout='c'
        h2.layout='c'
        q1.layout='c'
        b1.layout='c'
        b2.layout='c' 
        q2.layout='c'  
        w1.layout='c'
        w2.layout='c'         
        for m in range(self.m_min,self.m_max+1):
            m_local = m - self.m_min
            len_u1 = self.len_u1[m_local]
            len_h1 = self.len_h1[m_local]
            len_u2 = self.len_u2[m_local]
            len_h2 = self.len_h2[m_local]
            len_q1 = self.len_q1[m_local]
            len_b1 = self.len_b1[m_local]  
            len_b2 = self.len_b2[m_local]  
            len_q2 = self.len_q2[m_local]             
            len_w1 = self.len_w1[m_local]             
            u1['c'][m_local] = self.data[m_local][                                                 :len_u1]
            h1['c'][m_local] = self.data[m_local][len_u1                                           :len_u1+len_h1]
            u2['c'][m_local] = self.data[m_local][len_u1+len_h1                                    :len_u1+len_h1+len_u2]
            h2['c'][m_local] = self.data[m_local][len_u1+len_h1+len_u2                             :len_u1+len_h1+len_u2+len_h2]
            q1['c'][m_local] = self.data[m_local][len_u1+len_h1+len_u2+len_h2                      :len_u1+len_h1+len_u2+len_h2+len_q1]
            b1['c'][m_local] = self.data[m_local][len_u1+len_h1+len_u2+len_h2+len_q1               :len_u1+len_h1+len_u2+len_h2+len_q1+len_b1]
            b2['c'][m_local] = self.data[m_local][len_u1+len_h1+len_u2+len_h2+len_q1+len_b1        :len_u1+len_h1+len_u2+len_h2+len_q1+len_b1+len_b2]
            q2['c'][m_local] = self.data[m_local][len_u1+len_h1+len_u2+len_h2+len_q1+len_b1+len_b2 :len_u1+len_h1+len_u2+len_h2+len_q1+len_b1+len_b2+len_q2]
            w1['c'][m_local] = self.data[m_local][len_u1+len_h1+len_u2+len_h2+len_q1+len_b1+len_b2+len_q2 :len_u1+len_h1+len_u2+len_h2+len_q1+len_b1+len_b2+len_q2+len_w1]     
            w2['c'][m_local] = self.data[m_local][len_u1+len_h1+len_u2+len_h2+len_q1+len_b1+len_b2+len_q2+len_w1 :]                     
#def show_ball(S, field, index, longitude=0, latitude=0, mp = None):
#    
#    if mp == None:
#        figure, ax = plt.subplots(1,1)
#        figure.set_size_inches(3,3)
#
#    lon = np.linspace(0, 2*np.pi, 2*(S.L_max+1))
#    lat = S.grid - np.pi/2
#    
#    meshed_grid = np.meshgrid(lon, lat)
#    lat_grid    = meshed_grid[1]
#    lon_grid    = meshed_grid[0]
#    
#    if mp == None:
#        mp = Basemap(projection='ortho', lat_0=latitude, lon_0=longitude, ax=ax)
#        mp.drawmapboundary()
#        mp.drawmeridians(np.arange(0, 360, 30))
#        mp.drawparallels(np.arange(-90, 90, 30))
#
#    x, y = mp(np.degrees(lon_grid), np.degrees(lat_grid))
#    im = mp.pcolor(x, y, np.transpose(field), cmap='RdYlBu_r')
#    
#    
#    plt.savefig('images/om_%05i.png' %index)
#    return im,mp
