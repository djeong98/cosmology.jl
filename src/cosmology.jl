#!/usr/bin/env julia
#==============================================================================
>> Julia module for calculating cosmolgical quantities << 

 Original implementation in python (cosmology.py) 
 March 2013 by Donghui Jeong.
 Translated to Julia + Adding massive neutrino for the distance calculation
 For the equations and constants, see cosmo_eqns.pdf
 Oct. 2017 by Donghui Jeong

 Packaging 
 Jan. 2019 by Donghui Jeong

 =============================================================================#

#------------------------------------------------------------------------------
module cosmology # beginning of the module
#------------------------------------------------------------------------------

using Dierckx # For the spline interpolation
using OrdinaryDiffEq
using QuadGK			# integratin
using SpecialFunctions
using Printf
using DelimitedFiles

export cosParams,dAcz,dlnDdlna,Hz,Omz,Odez,Okz,Volume,Dplusz

const H0 = 0.00033356   			# in unit of [h/Mpc]
const speedOfLight = 3.0660e-7	# [Mpc/year]
const rho_crit = 2.7751e11			# [M_sun/h/(Mpc/h)^3]
const radian_in_degree = pi/180.
const radian_in_arcmin = radian_in_degree/60.
const radian_in_arcsec = radian_in_arcmin/60.
const Tcmb_pivot = 2.726			# in Kelvin
const Tnu_pivot = 0.16767e-3		# in eV
const Omegagh2 = 2.4746e-5			# photon energy density parameter
const OmegaNuMasslessPerNnuh2 = 5.7061e-6	# massless neutrino energy density parameter per d.o.f.
const OmegaNuMassiveWithoutfxnuh2 = 1.0042e-6	# massless neutrino energy density parameter per d.o.f.
const cm_in_pc = 3.0857e18			# parsec in cm

const xnumin = 1.e-4
const xnumax = 4000.
const fxnuInfCoeff = 1.5*zeta(3.)
const fnu0Val = 7/120*π^4.

struct cosParams{TP<:AbstractFloat}
	h::TP						# little h
	Om::TP						# Omega_m now
	Ocdm::TP					# Omeba_cdm now
	Ob::TP						# Omeba_b now
	Tnu::TP					# neutrino temperature now (in eV)
	aMnu::Array{TP,1} 	# neutrino mass vector (in eV)
	OnuMassive::TP			# Omega_nu now
	OnuMassless::TP			# Omega_nu now
	Ode::TP				# Omega_de now
	wde::TP					# w_de now
	Ok::TP					# Omega_curvature now
	Og::TP					# Omega_photo
	Hubble::TP			# prefactor for the distance (dimention 1/length)
	fxnu::Function	# neutrino energy density function
	fga::Spline1D		# growth of potential
	fdga::Spline1D		# derivative of ga
	fchiz::Spline1D		# chi(z)
	fzchi::Spline1D		# z(chi)
end

function cosParams{TP}(h=0.6778,Om=0.30821,Ode=0.69179;Ob=0.048555,wde=-1.,unit="Mpc/h",Tcmb=2.726,mnu=[0.06,0.,0.],calc_fxnu=false,calc_growth=false,calc_chiz=false) where {TP<:AbstractFloat}
	# Set the Unit
	if unit == "Mpc/h"
		Hubble = H0
	elseif unit == "Mpc"
		Hubble = H0*h
	else
		error("The unit must be either Mpc or Mpc/h")
	end

	# folder to store the pre-calculated data
	folder = "./tables"
	if isdir(folder)
		println("subdirectory "*folder*" exists!")
	else
		println("making the subdirectory "*folder*"!")
		mkdir(folder)
	end

	# For the massive neutrino energy density  
	fnufname = joinpath(folder,"frhonu_xnu.dat")
	if ~isfile(fnufname) || calc_fxnu
		calc_fxnu_function(fnufname)
	end
	fxnudata = readdlm(fnufname)
	fxnu_spline = Spline1D(fxnudata[:,1],fxnudata[:,2],k=5)
	fxnu = xnu -> (xnu>xnumax ? fxnuInfCoeff*xnu : (xnu<xnumin ? fnu0Val : fxnu_spline(xnu)))

	NnuMassless = count(mnu .== 0.)
	println("Number of massless neutrinos = $NnuMassless")
	NnuMassive  = 3-NnuMassless

	# Calculate the density parameters at present time
	# All radiation components are normalized with Tcmb. Set Tcmb=0 to turn them off.
	if Tcmb>0
		Og = Omegagh2*(Tcmb/Tcmb_pivot)^4/h^2
		Tnu = Tnu_pivot*(Tcmb/Tcmb_pivot)
		sumfxnu = sum([xnu>0 ? fxnu(xnu) : 0 for xnu in mnu/Tnu])
		OnuMassless = OmegaNuMasslessPerNnuh2*NnuMassless*(Tcmb/Tcmb_pivot)^4/h^2
		OnuMassive = OmegaNuMassiveWithoutfxnuh2*sumfxnu*(Tcmb/Tcmb_pivot)^4/h^2
		println("OnuMassless=$OnuMassless")
		println("OnuMassive=$OnuMassive")
	else
		Og = 0.
		OnuMassless  = 0.
	# If Tcmb is not given, using the default value for Tnu_pivot
		Tnu = Tnu_pivot
		sumfxnu = sum([xnu>0 ? fxnu(xnu) : 0 for xnu in mnu/Tnu])
		OnuMassive = OmegaNuMassiveWithoutfxnuh2*sumfxnu/h^2
	end

	Ocdm = Om - Ob - OnuMassive
	Or = Og + OnuMassless
	Ok = 1. - Om - Ode - Or

	println("Omegas at z=0 ==============================================")
	@printf "Omega_de           = %10.8f\n" Ode
	@printf "Omega_m            = %10.8f\n" Om
	@printf "Omega_cdm          = %10.8f\n" Ocdm
	@printf "Omega_b            = %10.8f\n" Ob
	@printf "Omega_nu(massive)  = %10.8f\n" OnuMassive
	@printf "Omega_g            = %10.8f\n" Og
	@printf "Omega_nu(massless) = %10.8f\n" OnuMassless
	@printf "Omega_k            = %10.8f\n" Ok
	println("============================================================")

	# the linear growth factor (Only using Omega_M, Omega_de, 1-Omega_M-Omdga_de at z=0)
	gfname = @sprintf("gfactor_Om=%6.4f_Ode=%6.4f_wde=%g.dat",Om,Ode,wde)
	gfname = joinpath(folder,gfname)
	if ~isfile(gfname) || calc_growth
		println("calling growth factor ")
		calc_growth_factor(gfname,Om,Ode,1-Om-Ode,wde)
	end
	gadata = readdlm(gfname)
	fga = Spline1D(gadata[:,1],gadata[:,2],k=5)
	fdga = Spline1D(gadata[:,1],gadata[:,3],k=5)

	# the chiz function
	if Tcmb>0
		chizfname = @sprintf("chiz_Om=%6.4f_Ode=%6.4f_wde=%g_Nmnu=%g_Mnu=%g_radiation_on.dat",Om,Ode,wde,NnuMassive,sum(mnu))
	else
		chizfname = @sprintf("chiz_Om=%6.4f_Ode=%6.4f_wde=%g_Nmnu=%g_Mnu=%g_radiation_off.dat",Om,Ode,wde,NnuMassive,sum(mnu))
	end
	chizfname = joinpath(folder,chizfname)
	if ~isfile(chizfname) | calc_chiz
		calc_chiz_function(chizfname,Ocdm+Ob,Ode,wde,Og+OnuMassless,Ok,fxnu,mnu,Tnu,1.0/Hubble,h)
	end
	chizdata =readdlm(chizfname)
	fchiz = Spline1D(chizdata[:,1],chizdata[:,2],k=5)
	fzchi = Spline1D(chizdata[:,2],chizdata[:,1],k=5)
	return cosParams(h,Om,Ocdm,Ob,Tnu,mnu,OnuMassive,OnuMassless,Ode,wde,Ok,Og,Hubble,fxnu,fga,fdga,fchiz,fzchi)
end
#==============================================================================
The functions needed for the constructors
==============================================================================#
#------------------------------------------------------------------------------
function calc_fxnu_function(nufname)
	println("Calculating the massive neutrino energy density function!")
	axnu = [10^x for x in range(log10(xnumin),stop=log10(xnumax),length=75)]
	afxnu = [quadgk(x-> x^2*sqrt(x^2+xnu^2)/(exp(x)+1),0,Inf,rtol=1.e-10)[1] for xnu in axnu]
	writedlm(nufname,zip(axnu,afxnu))
	return 0
end
#------------------------------------------------------------------------------
# t = lna; p=(Om,Ode,Ok,wde);u[1] = g; u[2] = dg/dlna
function fderivgrowth(du,u,p,t)
	a = exp(t)
	Hasq = p[1]/a^3 + p[2]/a^(3*(1.0+p[4])) + p[3]/a^2
	Oka  = p[3]/a^2/Hasq
	Odea = p[2]/a^(3*(1.0+p[4]))/Hasq
	du[1] = u[2]
	du[2] = -(2.5 + (Oka - 3.0*p[4]*Odea)/2.0)*u[2] - (2.0*Oka + 1.5*(1.0-p[4])*Odea)*u[1]
end
#------------------------------------------------------------------------------
function calc_growth_factor(gfname,Om,Ode,Ok,wde)

	alna = collect(range(-5,stop=0.01,length=502))

	println("Solving differential equation!")
#	fderivgrowth = PhiGrowth(Om=Om,Ode=Ode,Ok=Ok,wde=wde)

	u0    = [1.0;0.0]
	tspan = (-5.0,0.01)
    p = [Om,Ode,Ok,wde]
    prob  = ODEProblem(fderivgrowth,u0,tspan,p)
    @time sol = solve(prob,Tsit5(),saveat=alna,abstol=1.e-10,reltol=1.e-10)

	writedlm(gfname,zip(alna,sol[1,:],sol[2,:]))
	println("Done!!")

	return 0
end
#------------------------------------------------------------------------------
function calc_chiz_function(chizfname,Ocdmb,Ode,wde,Or,Ok,fxnu,mnu,Tnu,invH0,h)

	println("Integrating chiz")
	alogzp1 = collect(range(0,stop=3.5,length=351))
    function dchidx(x)
        fxnusum = sum([(xnu>0.0 ? fxnu(xnu) : 0.0) for xnu in mnu/(Tnu*x)])
        return 1.0/sqrt(Ocdmb*x^3 + Ode*x^(3(1+wde)) + Or*x^4 + Ok*x^2 + OmegaNuMassiveWithoutfxnuh2/h^2*x^4*fxnusum)
    end
    achiz   = [quadgk(dchidx,1,zp1,rtol=1.e-10)[1] for zp1 in 10 .^alogzp1]

	writedlm(chizfname,zip(10 .^alogzp1.-1,invH0.*achiz))
	println("Done!")
	return 0
end
#==============================================================================
The functions taking cosParams as an argument.
==============================================================================#
#------------------------------------------------------------------------------
function Rz(p::cosParams,z)
	return 3p.Ob/4p.Og/(1+z)
end
#------------------------------------------------------------------------------
function rhoMnu_ov_rhocrit(p::cosParams,a)
	return OmegaNuMassiveWithoutfxnuh2/p.h^2/a^4*sum([(xnu>0 ? p.fxnu(xnu) : 0) for xnu in p.aMnu/(p.Tnu/a)])
end
#------------------------------------------------------------------------------
function Hz(p::cosParams,z)
	zp1 = z+1
	return p.Hubble*sqrt( (p.Ocdm+p.Ob)*zp1^3 + p.Ode*zp1^(3(1+p.wde)) + (p.Og+p.OnuMassless)*zp1^4 + p.Ok*zp1^2 + rhoMnu_ov_rhocrit(p,1.0/(z+1)))
end
#------------------------------------------------------------------------------
function Ha(p::cosParams,a)
	return p.Hubble*sqrt( (p.Ocdm+p.Ob)/a^3 + p.Ode/a^(3(1+p.wde)) + (p.Og+p.OnuMassless)/a^4 + p.Ok/a^2 + rhoMnu_ov_rhocrit(p,a))
end
#------------------------------------------------------------------------------
function aHa(p::cosParams,a)
	return p.Hubble*sqrt( (p.Ocdm+p.Ob)/a + p.Ode/a^(1+3p.wde) + (p.Og+p.OnuMassless)/a^2 + p.Ok + rhoMnu_ov_rhocrit(p,a)*a^2 )
end
#------------------------------------------------------------------------------
function asqHa(p::cosParams,a)
	return p.Hubble*sqrt( (p.Ocdm+p.Ob)*a + p.Ode*a^(1-3p.wde) + (p.Og+p.OnuMassless) + p.Ok*a^2 + rhoMnu_ov_rhocrit(p,a)*a^4 )
end
#------------------------------------------------------------------------------
function Omz(p::cosParams,z)
	zp1 = z+1
	return ((p.Ocdm+p.Ob)*zp1^3 + rhoMnu_ov_rhocrit(p,1.0/zp1))/(Hz(p,z)/p.Hubble)^2
end
#------------------------------------------------------------------------------
function Odez(p::cosParams,z)
	zp1 = z+1
	return p.Ode*zp1^(3(1+p.wde))/(Hz(p,z)/p.Hubble)^2
end
#------------------------------------------------------------------------------
function Orz(p::cosParams,z)
	zp1 = z+1
	return (p.Og+p.OnuMassless)*zp1^4/(Hz(p,z)/p.Hubble)^2
end
#------------------------------------------------------------------------------
function Okz(p::cosParams,z)
	zp1 = z+1
	return p.Ok*zp1^2/(Hz(p,z)/p.Hubble)^2
end
#------------------------------------------------------------------------------
# Age of the Universe in year/h or year (depending on the unit)
#------------------------------------------------------------------------------
function agez(p::cosParams,z)
	a = 1.0/(z+1.0)
	agez = quadgk(x->1.0/aHa(p,x),0,a)[1]
	return agez/speedOfLight
end
#------------------------------------------------------------------------------
function comPtlHorizon(p::cosParams,z)
	a = 1.0/(z+1.0)
	return quadgk(x->1.0/asqHa(p,x),0,a)[1]
end
#------------------------------------------------------------------------------
function chiz(p::cosParams,z)
#	return quadgk(x->1./Hz(p,x),0,z)[1]
	return p.fchiz(z)
end
#------------------------------------------------------------------------------
function dAz(p::cosParams,z)
	if abs(p.Ok) < 1.e-4 
		return p.fchiz(z)/(1.0+z)
	elseif p.Ok > 0.0
		return sinh(p.fchiz(z)*sqrt(p.Ok)*p.Hubble)/(sqrt(p.Ok)*(1.0+z)*p.Hubble)
	elseif p.Ok < 0.0
		return sin(p.fchiz(z)*sqrt(-p.Ok)*p.Hubble)/(sqrt(-p.Ok)*(1.0+z)*p.Hubble)
	else
		error("If the code reaches here, something must be wrong. Check!")
	end
end
#------------------------------------------------------------------------------
function dAcz(p::cosParams,z)
	return (1.0+z)*dAz(p,z)
end
#------------------------------------------------------------------------------
function dLz(p::cosParams,z)
	return (1.0+z)^2*dAz(p,z)
end
#------------------------------------------------------------------------------
function Volume(p::cosParams,zi,zf,Area)
	return Area*radian_in_degree^2*quadgk(x->dAcz(p,x)^2/Hz(p,x),zi,zf)[1]
end
#------------------------------------------------------------------------------
function Dz(p::cosParams,z)
	lna = log(1.0/(1.0+z))
	ga = p.fga(lna)
	return ga/(1.0+z)
end
#------------------------------------------------------------------------------
function Dplusz(p::cosParams,z)
	return Dz(p,z)/Dz(p,0.0)
end
#------------------------------------------------------------------------------
function dlnDdlna(p::cosParams,z)
	lna = log(1.0/(1.0+z))
	ga = p.fga(lna)
	dga = p.fdga(lna)
	return 1.0+dga/ga
end
#------------------------------------------------------------------------------
function dfdlna(p::cosParams,z)
	a = 1.0/(1.0+z)
	Hasq = p.Om/a^3 + p.Ode/a^(3*(1.0+p.wde)) + (1.0-p.Om-p.Ode)/a^2
	Oka  = (1-p.Om-p.Ode)/a^2/Hasq
	Odea = p.Ode/a^(3*(1+p.wde))/Hasq
	lna = log(a)
	ga   = p.fga(lna)
	dga  = p.fdga(lna)
	ddga =  -(2.5 + (Oka - 3.0*p.wde*Odea)/2.0)*dga - (2.0*Oka + 1.5*(1.0-p.wde)*Odea)*ga
	return ddga/ga - (dlnDdlna(p,z)-1)^2
end
#------------------------------------------------------------------------------
function findz(p::cosParams,chi)
	return p.fzchi(chi)
end
#------------------------------------------------------------------------------
end # end of the module
#------------------------------------------------------------------------------
