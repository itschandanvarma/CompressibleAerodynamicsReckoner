import numpy as np
from scipy.optimize import fsolve

class Fluid:
 
    def __init__(self,gamma,R):
        self.gamma = gamma
        self.R = R

air = Fluid(1.4,286.9)

def SpeedofSound(Temperature,Fluid):
    a = np.sqrt(Fluid.gamma*Fluid.R*Temperature)
    return a

def MachNumber(Velocity,Temperature,Fluid):
    a = SpeedofSound(Temperature,Fluid)
    machnumber = Velocity/a
    return machnumber

#Stagnation Properties

class StagProp:

    def TstagTstatratio(MachNumber,Fluid):
        tsbyt = (1 + ((Fluid.gamma - 1)/2)*(MachNumber**2))
        return tsbyt

    def PstagPstatratio(MachNumber,Fluid):
        psbyp = StagProp.TstagTstatratio(MachNumber,Fluid)**(Fluid.gamma/(Fluid.gamma - 1))
        return psbyp

    def RhostagRhostatratio(MachNumber,Fluid):
        rhosbyrho = StagProp.TstagTstatratio(MachNumber,Fluid)**(1/(Fluid.gamma - 1))
        return rhosbyrho

    def Tstag(StaticTempearture,MachNumber,Fluid):
        ts = StaticTempearture*StagProp.TstagTstatratio(MachNumber,Fluid)
        return ts

    def Pstag(StaticPressure,MachNumber,Fluid):
        ps = StaticPressure*StagProp.PstagPstatratio(MachNumber,Fluid)
        return ps

    def Rhostag(StaticDensity,MachNumber,Fluid):
        rhos = StaticDensity*StagProp.RhostagRhostatratio(MachNumber,Fluid)
        return rhos

#Sonic Properties

class SonicProp:

    def TstagTstarratio(Fluid):
        MachNumber = 1
        tsbytstar = (1 + ((Fluid.gamma - 1)/2)*(MachNumber**2))
        return tsbytstar

    def PstagPstarratio(Fluid):
        MachNumber = 1
        psbypstar = StagProp.TstagTstatratio(MachNumber,Fluid)**(Fluid.gamma/(Fluid.gamma - 1))
        return psbypstar

    def RhostagRhostarratio(Fluid):
        MachNumber = 1
        rhosbyrhostar = StagProp.TstagTstatratio(MachNumber,Fluid)**(1/(Fluid.gamma - 1))
        return rhosbyrhostar

    def TstatTstarratio(MachNumber,Fluid):
        tbytstar = SonicProp.TstagTstarratio(Fluid)/StagProp.TstagTstatratio(MachNumber,Fluid)
        return tbytstar

    def PstatPstarratio(MachNumber,Fluid):
        pbypstar = SonicProp.PstagPstarratio(Fluid)/StagProp.PstagPstatratio(MachNumber,Fluid)
        return pbypstar

    def RhostatRhostarratio(MachNumber,Fluid):
        rhobyrhostar = SonicProp.RhostagRhostarratio(Fluid)/StagProp.RhostagRhostatratio(MachNumber,Fluid)
        return rhobyrhostar

    def AstatAstarratio(MachNumber,Fluid):
        abyastar = np.sqrt(SonicProp.TstagTstarratio(Fluid)/StagProp.TstagTstatratio(MachNumber,Fluid))
        return abyastar

    def Mstar(MachNumber,Fluid):
        mstar = MachNumber*(np.sqrt(SonicProp.TstagTstarratio(Fluid)/StagProp.TstagTstatratio(MachNumber,Fluid)))
        return mstar

#Normal Shock Relations

class NormalShocks:

    def Rho2byRho1(M1,Fluid):
        rho2byrho1 = np.square(SonicProp.Mstar(M1,Fluid))
        return rho2byrho1

    def P2byP1(M1,Fluid):
        p2byp1 = 1 + (Fluid.gamma*(M1**2)*(1 - (1/NormalShocks.Rho2byRho1(M1,Fluid))))
        return p2byp1

    def T2byT1(M1,Fluid):
        t2byt1 = NormalShocks.P2byP1(M1,Fluid)/NormalShocks.Rho2byRho1(M1,Fluid)
        return t2byt1

    def M2(M1,Fluid):
        m2 = np.sqrt(StagProp.TstagTstatratio(M1,Fluid)/((Fluid.gamma*(M1**2) - (0.5*(Fluid.gamma - 1)))))
        return m2

    def P02byP01(M1,Fluid):
        M2 = NormalShocks.M2(M1,Fluid)
        p02byp01 = (StagProp.PstagPstatratio(M2,Fluid)/StagProp.PstagPstatratio(M1,Fluid))*(NormalShocks.P2byP1(M1,Fluid))
        return p02byp01

#Oblique Shock Relations

class ObliqueShocks:
        
    def DeflectionAngle(M,ShockAngle,Fluid):
        beta = np.deg2rad(ShockAngle)
        theta = np.arctan((2*(1/np.tan(beta)))*((((M*np.sin(beta))**2) - 1)/(((M**2)*(Fluid.gamma + np.cos(2*beta))) + 2)))
        theta = np.rad2deg(theta)
        return theta
    
    def ShockAngle(M,DeflectionAngle,Fluid):
        a = (1 + (0.5*(Fluid.gamma - 1)*(M**2)))*np.tan(np.deg2rad(DeflectionAngle))
        b = 1 - (M**2)
        c = (1 + (0.5*(Fluid.gamma + 1)*(M**2)))*np.tan(np.deg2rad(DeflectionAngle))
        d = 1
        roots = np.roots([a,b,c,d])
        roots = np.rad2deg(np.arctan(roots))
        WeakShockAngle = min(roots[0],roots[1])
        StrongShockAngle = max(roots[0],roots[1])
        return WeakShockAngle,StrongShockAngle

    def Rho2byRho1(M1,ShockAngle,Fluid):
        beta = np.deg2rad(ShockAngle)
        mn1 = M1*np.sin(beta)
        rho2byrho1 = np.square(SonicProp.Mstar(mn1,Fluid))
        return rho2byrho1

    def P2byP1(M1,ShockAngle,Fluid):
        beta = np.deg2rad(ShockAngle)
        mn1 = M1*np.sin(beta)
        p2byp1 = 1 + (Fluid.gamma*(mn1**2)*(1 - (1/NormalShocks.Rho2byRho1(mn1,Fluid))))
        return p2byp1

    def T2byT1(M1,ShockAngle,Fluid):
        beta = np.deg2rad(ShockAngle)
        mn1 = M1*np.sin(beta)
        t2byt1 = NormalShocks.P2byP1(mn1,Fluid)/NormalShocks.Rho2byRho1(mn1,Fluid)
        return t2byt1

    def M2(M1,ShockAngle,Fluid):
        beta = np.deg2rad(ShockAngle)
        mn1 = M1*np.sin(beta)
        mn2 = np.sqrt(StagProp.TstagTstatratio(mn1,Fluid)/((Fluid.gamma*(mn1**2) - (0.5*(Fluid.gamma - 1)))))
        theta = ObliqueShocks.DeflectionAngle(M1,ShockAngle,air)
        theta = np.deg2rad(theta)
        m2 = mn2/np.sin(beta-theta)
        return m2

    def P02byP01(M1,ShockAngle,Fluid):
        beta = np.deg2rad(ShockAngle)
        mn1 = M1*np.sin(beta)
        M2 = ObliqueShocks.M2(M1,ShockAngle,Fluid)
        p02byp01 = (StagProp.PstagPstatratio(M2,Fluid)/StagProp.PstagPstatratio(mn1,Fluid))*(NormalShocks.P2byP1(mn1,Fluid))
        return p02byp01

#Prandtl Meyer Relations

class PrandtlMeyer:

    def PrandtlMeyerAngle(M,Fluid):
        nu = np.sqrt((Fluid.gamma + 1)/(Fluid.gamma - 1))*np.rad2deg(np.arctan(np.sqrt((Fluid.gamma - 1)*((M**2) - 1)/(Fluid.gamma + 1)))) - np.rad2deg(np.arctan(np.sqrt((M**2) - 1)))
        return nu

    def MachNumber(PrandtlMeyerAngle,Fluid):
        func = lambda M : np.sqrt((Fluid.gamma + 1)/(Fluid.gamma - 1))*np.rad2deg(np.arctan(np.sqrt((Fluid.gamma - 1)*((M**2) - 1)/(Fluid.gamma + 1)))) - np.rad2deg(np.arctan(np.sqrt((M**2) - 1))) - PrandtlMeyerAngle
        M = fsolve(func,1)
        return M