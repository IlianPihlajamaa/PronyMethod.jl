using Test, Prony


@testset "Prony" begin
    for N = [30, 31]
        x = collect(LinRange(0, 10, N))
        y = @. exp(-x*2)*cos(3x)+sin(2x)*exp(-x*3)
        f  = prony(x,y)
        ynew = real.(f.(x))
        error = maximum(abs.(y.-ynew))
        @test error < 10^-1
    end
end

@testset "ApproximateProny" begin
    for N = [30, 31]
        x = collect(LinRange(0, 4, N))
        y = @. exp(-x*2)*cos(3x)+sin(2x)*exp(-x*3)
        f  = prony(x,y, 10)
        ynew = real.(f.(x))
        error = maximum(abs.(y.-ynew))
        @test error < 10^-1
    end
end