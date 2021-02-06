using SafeTestsets
using Test


@safetestset vas_sample "VAS samples with direct" begin
    rng = MersenneTwister(2930472)
    vas = VectorAdditionModel(transitions, rates)
    sampler = MarkovDirect()

    input_process = vas_initial(vas, [1, 1, 0])
    for i in 1:10
        y = send(vas, input_process)
        sample = send(sampler, y)
        input_process = vas_input(vas, sample)
    end
end
