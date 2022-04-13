using Test, Prony


@testset begin
    x = collect(LinRange(0, 5, 50))
    y = @. exp(-x*2)*cos(3x)+sin(2x)
    for f in [prony(x,y), prony(x,y, 20), prony(x,y, 20, PronyMethodMPM())]
        ynew = f.(x)
        error = maximum(abs.(y.-ynew))
        @test error < 10^-5
    end
end