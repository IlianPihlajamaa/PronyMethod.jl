using Test, Prony


@testset "Prony" begin
    for N = 20:10:60
        x = collect(LinRange(0, 4, N))
        y = @. exp(-x*2)*cos(3x)+sin(2x)
        f  = prony(x,y)
        ynew = real.(f.(x))
        error = maximum(abs.(y.-ynew))
        @test error < 10^-1
    end
end

@testset "ApproximateProny" begin
    for N = 20:10:60
        x = collect(LinRange(0, 4, N))
        y = @. exp(-x*2)*cos(3x)+sin(2x)*exp(-x*3)
        f  = prony(x,y, 10)
        ynew = real.(f.(x))
        error = maximum(abs.(y.-ynew))
        @test error < 10^-1
    end
end