using NoisyNeuron, Test

println("Let's test things")

@testset "point" begin
  @test length(point_white_steady_test()) == 2
  @test length(point_white_cmod_test()) == 2
  @test length(point_white_vmod_test()) == 2
end

@testset "cable" begin
  @test length(sealed_white_steady_var_test()) == 3
  @test length(sealed_colour_steady_var_test()) == 5
end