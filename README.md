# Bayesian VaR Backtester

An interactive app for learning Bayesian updating through Value-at-Risk backtesting. It tracks a 99% one-day VaR model's breach rate over a 250-day backtest year. After each 25-day window, it forecasts the year-end Basel traffic-light zone and shows how uncertain that forecast is.

The app is a teaching tool, not a validation or capital-planning system. It favors clear, visible math over regulatory precision.

**Live app:** `https://patrick-lefler.com/bayesian_var_backtester/`

## What the App Shows

After each 25-day window, the Simulator tab shows four results:

1. **25-day breach rate**: breaches divided by trading days in that window.
2. **Year-to-date rate**: all breaches and trading days so far, plus the current zone.
3. **Bayes forecast of the Basel zone at day 250**: the probability of ending green, yellow, or red, and the expected breach count with a 95% range.
4. **Uncertainty curve**: the full distribution of year-end breaches, colored by zone and layered over the curves from earlier windows.

It also shows the posterior estimate of the true breach rate θ and the probability that θ is above 1%. Other panels include a clustering check, a step-by-step breakdown of the calculation, a zone-forecast chart, a breach-rate trajectory, and an editable backtest log.

The Guide & Math tab walks through the model. Its worked example updates live from whatever values you enter.

## The Model

| Piece | Specification |
|---|---|
| True breach rate | θ ~ Beta(α₀, β₀) |
| Observed breaches | B \| θ ~ Binomial(D, θ) |
| Updated belief | θ ~ Beta(α₀ + B, β₀ + D − B) |
| Remaining breaches | Beta-Binomial(250 − D, α, β) |
| Year-end count | B + remaining breaches |

**Prior.** The prior comes from the two previous backtest years, each discounted by a weight:

- α₀ = w₁·B₁ + w₂·B₂
- β₀ = w₁·(D₁ − B₁) + w₂·(D₂ − B₂)

By default, last year counts for 0.5 and two years ago for 0.4. With 3 and 2 breaches in those years, the prior is Beta(2.2, 222.8): a mean of 0.98%, worth 225 trading days of evidence. Presets let you compare it with pooling both years, using last year only, or the "Model's claim" prior. That prior is centered on exactly 1%, and a slider sets how many days of evidence it's worth (default 500).

**Zones.** The app uses the 1996 Basel traffic light for 250 days at 99% VaR:

| Zone | Breaches | Multiplier |
|---|---|---|
| Green | 0–4 | 3.00 |
| Yellow | 5–9 | 3.40, 3.50, 3.65, 3.75, 3.85 |
| Red | 10 or more | 4.00 |

**Uncertainty.** The 95% range uses exact Beta-Binomial quantiles. The guide splits the variance of remaining breaches into two parts:

- **Luck:** R·p·(1 − p), the randomness of the days still to come.
- **Rate uncertainty:** R²·Var(θ), the remaining doubt about the model's true breach rate. This part shrinks as the year goes on.

**Clustering caveat.** The Beta-Binomial forecast assumes breaches are independent. A rough dispersion test compares the 25-day window counts with what independence would produce. The result, D = Σ(x − 25p̂)² / (25p̂(1 − p̂)), is checked against χ² with k − 1 degrees of freedom. The app flags "Possible clustering" at p < 0.10 and "Clustering detected" at p < 0.05, once there are at least 3 windows and 3 breaches. With so few breaches the test has little power. It's a prompt to run a proper independence test, such as Christoffersen's, not a substitute for one.

## Using the App

- Edit the breaches and trading days for either prior year, or open "Prior source & weights" to change the weights or pick a preset.
- Step through the year with the arrows or press Play. Clicking a row in the backtest log jumps to that window.
- Edit any 25-day window to see how a cluster of breaches moves the zone forecast and the clustering flag.

**Default data.** The current year's windows are [0, 1, 0, 0, 0, 3, 3, 0, 0, 0]. The model sits in green through day 150, then six breaches arrive in two windows. It ends the year with 7 breaches, in yellow, and the clustering flag trips.

**Input rules:**

- Breaches and trading days must be whole numbers, with 0 ≤ breaches ≤ days. Each 25-day window allows 0–25 breaches.
- Each prior year needs at least 1 trading day unless its weight is 0.
- α₀ and β₀ are floored at 0.5 so the Beta distribution stays valid.
- Prior years above 260 trading days get a soft warning.
- More than 5 breaches in one 25-day window gets a soft warning.

## Running Locally

There's no build step and nothing to install. Clone the repo and open `index.html` in a browser.

```
index.html   page structure and static guide content
styles.css   dark theme and layout
app.js       Beta-Binomial math, zone forecast, dispersion test, SVG charts, and dynamic guide text
preview.png  link-preview image
```

The app is plain JavaScript with hand-drawn SVG charts and no third-party libraries. Fonts load from Google Fonts. Without a connection, the browser falls back to system fonts.

## Sources and Credits

**Regulatory framework**

- Basel Committee on Banking Supervision (1996), ["Supervisory framework for the use of 'backtesting' in conjunction with the internal models approach to market risk capital requirements"](https://www.bis.org/publications/199601-standards-supervisory-framework-use-backtesting-conjunction-internal-models-approach-market-risk-capital.pdf). It's the source of the green, yellow, and red zones, the 250-day window, and the plus factors.
- Basel Framework, [MAR32](https://www.bis.org/committees/bcbs/basel-framework/81326/chapter.pdf) and [MAR99](https://www.bis.org/committees/bcbs/basel-framework/81307/chapter.pdf). The current standard keeps the 250-day, 99% count but calls the middle zone "amber" and uses a different multiplier scale under the Fundamental Review of the Trading Book.

**Backtesting methods**

- Paul H. Kupiec (1995), ["Techniques for Verifying the Accuracy of Risk Measurement Models,"](https://doi.org/10.3905/jod.1995.407942) Journal of Derivatives 3(2): 73–84. See also the [Federal Reserve working paper](https://fraser.stlouisfed.org/title/finance-economics-discussion-series-1491/techniques-verifying-accuracy-risk-measurement-models-717921).
- Peter F. Christoffersen (1998), ["Evaluating Interval Forecasts,"](https://doi.org/10.2307/2527341) International Economic Review 39(4): 841–862.

**Companion app**

- [Bayes at the Plate](https://patrick-lefler.com/bayes_batting_avg_learning_app/) uses the same Beta-Binomial engine, applied to baseball batting averages.

**Fonts:** [Inter](https://fonts.google.com/specimen/Inter) and [JetBrains Mono](https://fonts.google.com/specimen/JetBrains+Mono), used under the SIL Open Font License.

## Disclaimer

- The default backtest data is illustrative and doesn't represent any real bank, desk, or model.
- This project isn't affiliated with or endorsed by the Basel Committee, the BIS, any supervisor, or any author cited above.
- It's an educational tool, not regulatory, compliance, or capital-planning advice.

## Author

Patrick Lefler · Seeking Eureka

## License

MIT. See `LICENSE.md`.
