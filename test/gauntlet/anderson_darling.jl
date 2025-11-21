using HypothesisTests


"""
  - High p-value (> 0.05): Fail to reject H₀ → samples are consistent with coming from the same
  distribution ✓ (good for calibration)
  - Low p-value (≤ 0.05): Reject H₀ → evidence that samples come from different distributions ✗
  (indicates a problem)
"""
function ad_two_sample(draws_a::Vector{ClockDraw}, draws_b::Vector{ClockDraw}, clocks)
    # Test holding times which are times for specific clocks to fire.
    for clock in clocks
        times_a = [x[2] for x in draws_a if x[1] == clock]
        times_b = [x[2] for x in draws_b if x[1] == clock]
        result = KSampleADTest(times_a, times_b)
        @show clock, pvalue(result)
    end
    # Test waiting time which is over all clocks.
    times_a = [x[2] for x in draws_a]
    times_b = [x[2] for x in draws_b]
    result = KSampleADTest(times_a, times_b)
    @show pvalue(result)
end
