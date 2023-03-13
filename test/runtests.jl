using gravitywaveslagrangian
using Test

@testset "gravitywaveslagrangian.jl" begin
    @test gravitywaveslagrangian.test_cosmology() == 2.3e-10
    println("Cosmology test passed.")
end
