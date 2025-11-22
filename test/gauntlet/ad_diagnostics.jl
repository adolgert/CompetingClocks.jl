"""
    ad_diagnostic_report(test::KSampleADTest)

Provides a detailed diagnostic report for Anderson-Darling test results,
helping assess the reliability and interpretation of the test.
"""
function ad_diagnostic_report(test::HypothesisTests.KSampleADTest)
    println("\n" * "="^60)
    println("Anderson-Darling Test Diagnostic Report")
    println("="^60)

    # 1. Sample sizes
    println("\n1. SAMPLE SIZES:")
    for (i, size) in enumerate(test.sizes)
        println("   Sample $i: $size observations")
    end
    total = sum(test.sizes)
    println("   Total: $total observations")

    # Sample size adequacy
    min_size = minimum(test.sizes)
    if min_size < 30
        println("   ⚠️  WARNING: Small sample size (n < 30). Test may lack power.")
        println("      Consider using nsim parameter for simulated p-values.")
    elseif min_size < 100
        println("   ⚡ Moderate sample size. Test has reasonable power.")
    else
        println("   ✓ Large sample size. Test has good power.")
    end

    # 2. Test statistic
    println("\n2. TEST STATISTIC:")
    println("   A²k = $(test.A²k)")
    println("   SD(A²k) = $(test.σ)")
    println("   Standardized = $(test.A²k / test.σ)")

    # Interpretation of statistic magnitude
    std_stat = test.A²k / test.σ
    if std_stat < 1.0
        println("   → Small deviation from H₀ (< 1 SD)")
    elseif std_stat < 2.0
        println("   → Moderate deviation from H₀ (1-2 SD)")
    else
        println("   → Large deviation from H₀ (> 2 SD)")
    end

    # 3. P-value and interpretation
    pval = pvalue(test)
    println("\n3. P-VALUE:")
    println("   p-value = $pval")
    println("   Method: $(test.nsim == 0 ? "asymptotic" : "simulated (nsim=$(test.nsim))")")

    if pval < 0.001
        println("   *** Strong evidence of difference (p < 0.001)")
    elseif pval < 0.01
        println("   **  Strong evidence of difference (p < 0.01)")
    elseif pval < 0.05
        println("   *   Evidence of difference (p < 0.05)")
    elseif pval < 0.10
        println("   ⚡  Weak evidence of difference (0.05 < p < 0.10)")
    else
        println("   ✓   No evidence of difference (p ≥ 0.10)")
    end

    # 4. Recommendations
    println("\n4. RECOMMENDATIONS:")

    if min_size < 30 && test.nsim == 0
        println("   • Use nsim=1000 or higher for more reliable p-values with small samples")
        println("     Example: KSampleADTest(s1, s2, nsim=1000)")
    end

    if test.nsim > 0 && test.nsim < 1000
        println("   • Consider nsim ≥ 1000 for stable p-value estimates")
    end

    # Check if we're doing multiple tests
    println("   • If testing multiple clocks, apply multiple testing correction")
    println("     - Bonferroni: use α/k where k = number of tests")
    println("     - For 5 tests at α=0.05, use threshold p < 0.01")

    # Power considerations
    if pval > 0.05 && min_size < 100
        println("   • Non-significant result with moderate sample size:")
        println("     - May lack power to detect small differences")
        println("     - Consider increasing sample size if possible")
        println("     - Check effect size (A²k statistic magnitude)")
    end

    println("\n" * "="^60)
end


"""
    ad_summary_table(results::Vector{Tuple{Any, KSampleADTest}})

Creates a summary table for multiple Anderson-Darling tests,
useful when testing multiple clocks or conditions.
"""
function ad_summary_table(results::Vector)
    println("\n" * "="^80)
    println("Anderson-Darling Test Summary")
    println("="^80)
    println(rpad("Test", 15), rpad("n₁", 8), rpad("n₂", 8),
            rpad("A²k", 12), rpad("p-value", 12), "Significant?")
    println("-"^80)

    significant_count = 0
    bonferroni_alpha = 0.05 / length(results)

    for (label, test) in results
        pval = pvalue(test)
        sig = pval < 0.05 ? "*" : ""
        bonf_sig = pval < bonferroni_alpha ? "**" : sig

        println(rpad(string(label), 15),
                rpad(string(test.sizes[1]), 8),
                rpad(string(test.sizes[2]), 8),
                rpad(string(round(test.A²k, digits=4)), 12),
                rpad(string(round(pval, digits=4)), 12),
                bonf_sig)

        if pval < bonferroni_alpha
            significant_count += 1
        end
    end

    println("-"^80)
    println("* p < 0.05 (uncorrected)")
    println("** p < $(round(bonferroni_alpha, digits=4)) (Bonferroni corrected for $(length(results)) tests)")
    println("\nSignificant tests (Bonferroni): $significant_count / $(length(results))")
    println("="^80)
end
