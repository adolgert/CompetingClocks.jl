using HypothesisTests


"""
  - High p-value (> 0.05): Fail to reject H₀ → samples are consistent with coming from the same
  distribution ✓ (good for calibration)
  - Low p-value (≤ 0.05): Reject H₀ → evidence that samples come from different distributions ✗
  (indicates a problem)

Set verbose=true to get detailed diagnostic information about test quality and reliability.
"""
function ad_two_sample(draws_a::Vector{ClockDraw}, draws_b::Vector{ClockDraw}, clocks; verbose=false)
    results = Any[]

    # Test holding times which are times for specific clocks to fire.
    for clock in clocks
        times_a = [x[2] for x in draws_a if x[1] == clock]
        times_b = [x[2] for x in draws_b if x[1] == clock]
        result = KSampleADTest(times_a, times_b)
        pv = pvalue(result)
        push!(results, (;test="ad-two-sample", clock, pvalue=pv, result))
        @show clock, pv

        if verbose
            ad_diagnostic_report(result)
        end
    end

    # Test waiting time which is over all clocks.
    times_a = [x[2] for x in draws_a]
    times_b = [x[2] for x in draws_b]
    result = KSampleADTest(times_a, times_b)
    pv = pvalue(result)
    push!(results, (; test="ad-two-sample", clock=0, pvalue, result))

    if verbose
        ad_diagnostic_report(result)
        ad_summary_table(results)
    end

    return results
end
