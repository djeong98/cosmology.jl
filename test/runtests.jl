using Test
using cosmology

@testset "cosmology.jl" begin

    @testset "Initialization - Planck 2018" begin
        # Test that we can initialize with Planck 2018 parameters
        cosmo = cosParams(
            h = 0.6766,
            Om = 0.3111,
            Ode = 0.6889,
            Ob = 0.048975,
            wde = -1.0,
            unit = "Mpc/h",
            Tcmb = 2.726
        )

        # Check that basic parameters are stored correctly
        @test cosmo.h ≈ 0.6766
        @test cosmo.Om ≈ 0.3111
        @test cosmo.Ode ≈ 0.6889
        @test cosmo.Ob ≈ 0.048975
        @test cosmo.wde ≈ -1.0

        # Check that density parameters sum to ~1
        total_Omega = cosmo.Om + cosmo.Ode + cosmo.Ok + cosmo.Og + cosmo.OnuMassless
        @test total_Omega ≈ 1.0 atol=1e-6

        println("✓ Initialization tests passed")
    end

    @testset "Hubble Parameter" begin
        cosmo = cosParams(h=0.7, Om=0.3, Ode=0.7, Tcmb=0.0)

        # At z=0, H(z) should equal H0
        H0_expected = 0.00033356 * 0.7  # in 1/(Mpc/h)
        @test Hz(cosmo, 0.0) ≈ H0_expected rtol=1e-10

        # H(z) should increase with redshift
        H_z1 = Hz(cosmo, 1.0)
        H_z2 = Hz(cosmo, 2.0)
        @test H_z2 > H_z1 > H0_expected

        # Test consistency between Hz and Ha
        z = 1.5
        a = 1.0 / (1.0 + z)
        @test Hz(cosmo, z) ≈ Ha(cosmo, a) rtol=1e-10

        println("✓ Hubble parameter tests passed")
    end

    @testset "Density Parameters Evolution" begin
        cosmo = cosParams(h=0.7, Om=0.3, Ode=0.7, Tcmb=0.0)

        # At z=0
        @test Omz(cosmo, 0.0) ≈ 0.3 rtol=1e-10
        @test Odez(cosmo, 0.0) ≈ 0.7 rtol=1e-10
        @test Okz(cosmo, 0.0) ≈ 0.0 atol=1e-10

        # Matter density should increase with redshift
        @test Omz(cosmo, 1.0) > Omz(cosmo, 0.0)
        @test Omz(cosmo, 2.0) > Omz(cosmo, 1.0)

        # Dark energy density should decrease with redshift (for w=-1)
        @test Odez(cosmo, 1.0) < Odez(cosmo, 0.0)

        # All density parameters should sum to 1
        z = 1.5
        total = Omz(cosmo, z) + Odez(cosmo, z) + Okz(cosmo, z) + Orz(cosmo, z)
        @test total ≈ 1.0 atol=1e-6

        println("✓ Density parameter evolution tests passed")
    end

    @testset "Distance Functions" begin
        cosmo = cosParams(h=0.7, Om=0.3, Ode=0.7, Tcmb=0.0)

        # Comoving distance should be zero at z=0
        @test chiz(cosmo, 0.0) ≈ 0.0 atol=1e-10

        # Comoving distance should increase with redshift
        chi1 = chiz(cosmo, 1.0)
        chi2 = chiz(cosmo, 2.0)
        @test chi2 > chi1 > 0

        # Check distance relations
        z = 1.5
        dA = dAz(cosmo, z)
        dL = dLz(cosmo, z)
        dAc = dAcz(cosmo, z)

        # dAc = (1+z) * dA
        @test dAc ≈ (1+z) * dA rtol=1e-10

        # dL = (1+z)^2 * dA
        @test dL ≈ (1+z)^2 * dA rtol=1e-10

        # For flat universe, dAc ≈ chi
        @test dAc ≈ chiz(cosmo, z) rtol=1e-6

        # Test findz - should be inverse of chiz
        chi_test = chiz(cosmo, 1.0)
        z_recovered = findz(cosmo, chi_test)
        @test z_recovered ≈ 1.0 rtol=1e-6

        println("✓ Distance function tests passed")
    end

    @testset "Growth Functions" begin
        cosmo = cosParams(h=0.7, Om=0.3, Ode=0.7, Tcmb=0.0)

        # Normalized growth should be 1 at z=0
        @test Dplusz(cosmo, 0.0) ≈ 1.0 rtol=1e-10

        # Growth factor should decrease with increasing redshift
        D_z1 = Dplusz(cosmo, 1.0)
        D_z2 = Dplusz(cosmo, 2.0)
        @test D_z2 < D_z1 < 1.0

        # dlnD/dlna should be positive
        f = dlnDdlna(cosmo, 1.0)
        @test f > 0

        # At high z (matter dominated), f ≈ 1
        f_highz = dlnDdlna(cosmo, 10.0)
        @test f_highz ≈ 1.0 rtol=0.1  # Within 10%

        println("✓ Growth function tests passed")
    end

    @testset "Volume Calculation" begin
        cosmo = cosParams(h=0.7, Om=0.3, Ode=0.7, Tcmb=0.0)

        # Volume should be positive
        V = Volume(cosmo, 0.5, 1.5, 1.0)  # 1 deg²
        @test V > 0

        # Volume should increase with survey area
        V1 = Volume(cosmo, 1.0, 2.0, 10.0)
        V2 = Volume(cosmo, 1.0, 2.0, 20.0)
        @test V2 ≈ 2 * V1 rtol=1e-10

        # Volume should increase with redshift range
        V_small = Volume(cosmo, 1.0, 1.5, 1.0)
        V_large = Volume(cosmo, 1.0, 2.0, 1.0)
        @test V_large > V_small

        println("✓ Volume calculation tests passed")
    end

    @testset "Unit System" begin
        # Test Mpc/h vs Mpc units
        cosmo_h = cosParams(h=0.7, Om=0.3, Ode=0.7, unit="Mpc/h", Tcmb=0.0)
        cosmo_nonh = cosParams(h=0.7, Om=0.3, Ode=0.7, unit="Mpc", Tcmb=0.0)

        z = 1.0

        # Distances in Mpc should be h times larger than in Mpc/h
        chi_h = chiz(cosmo_h, z)
        chi_nonh = chiz(cosmo_nonh, z)
        @test chi_nonh ≈ chi_h / 0.7 rtol=1e-6

        println("✓ Unit system tests passed")
    end

    @testset "Massive Neutrinos" begin
        # Cosmology without massive neutrinos
        cosmo_nomnu = cosParams(h=0.7, Om=0.3, Ode=0.7, mnu=[0.0, 0.0, 0.0], Tcmb=0.0)

        # Cosmology with one massive neutrino
        cosmo_mnu = cosParams(h=0.7, Om=0.3, Ode=0.7, mnu=[0.06, 0.0, 0.0], Tcmb=0.0)

        # Massive neutrino should increase matter density
        @test cosmo_mnu.OnuMassive > 0
        @test cosmo_nomnu.OnuMassive ≈ 0 atol=1e-10

        # This affects distances at high redshift
        z_high = 3.0
        chi_nomnu = chiz(cosmo_nomnu, z_high)
        chi_mnu = chiz(cosmo_mnu, z_high)
        @test abs(chi_mnu - chi_nomnu) / chi_nomnu > 0.001  # >0.1% difference

        println("✓ Massive neutrino tests passed")
    end

    @testset "Extreme Cases" begin
        cosmo = cosParams(h=0.7, Om=0.3, Ode=0.7, Tcmb=0.0)

        # Very small redshift
        z_small = 0.001
        @test Hz(cosmo, z_small) > 0
        @test chiz(cosmo, z_small) > 0

        # High redshift
        z_high = 10.0
        @test Hz(cosmo, z_high) > 0
        @test chiz(cosmo, z_high) > 0

        println("✓ Extreme case tests passed")
    end

end

println("\n" * "="^60)
println("All cosmology.jl tests completed successfully!")
println("="^60)
