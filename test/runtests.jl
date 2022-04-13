using Test, Prony, SpecialFunctions


@testset "Prony" begin
    for N = 50:50:200
        x = collect(LinRange(0, 4, 50))
        y = @. besselj0(x)#exp(-x*2)*cos(3x)+sin(2x)
        f  = prony(x,y)
        ynew = f.(x)
        error = maximum(abs.(y.-ynew))
        @test error < 10^-5
    end
end

@testset "PronyLS" begin
    for N = 100:50:200
        x = collect(LinRange(0, 5, N))
        y = @. besselj0(x)
        f  = prony(x,y,round(Int64,N/2-2))
        ynew = f.(x)
        error = maximum(abs.(y.-ynew))
        @test error < 10^-5
    end
end

@testset "PronyMPM" begin
    for N = 100:50:200
        x = collect(LinRange(0, 5, N))
        y = @. besselj0(x)
        f  = prony(x,y,round(Int64,N/2-10), PronyMethodMPM())
        ynew = f.(x)
        error = maximum(abs.(y.-ynew))
        @test error < 10^-5
    end
end