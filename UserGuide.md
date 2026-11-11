# Bayesian VaR Backtester: A Plain-English Guide

## The Big Idea

Banks use a risk model called Value at Risk, or VaR, to answer a simple question: "How much could we lose on a bad day?" A 99% VaR means that on 99 out of 100 trading days, losses should stay below a certain amount.

On the remaining 1 day out of 100, losses go above that amount. This is called a breach. Breaches are normal. In a year with 250 trading days, a good model should have about 2 or 3 breaches.

Regulators track breaches using a traffic light system:

| Zone | Breaches in 250 days | What it means |
|---|---|---|
| Green | 0 to 4 | The model looks fine. |
| Yellow | 5 to 9 | Questions get asked, and the bank may have to hold more capital. |
| Red | 10 or more | The model is presumed broken, and the capital penalty is largest. |

The issue is that the score only matters at year-end. If you see 4 breaches halfway through, it still looks “green,” but that doesn’t really show where things will end up. This app helps answer what boards really want to know: **“Where are we likely to finish, and how confident can we be?”**

## How the App Thinks

The app works in a practical way. It begins with what it already knows and updates its view as new information comes in. This approach is called Bayes’ Theorem.

1. **Start with history.** The model had 3 breaches two years ago and 2 last year. Last year counts a little more because it’s more recent.
2. **Add new evidence.** Every 25 trading days, or about a month, the app updates the record with the latest breaches.
3. **Forecast the finish.** The app figures out the chances of ending the year in green, yellow, or red.

At the start of the year, history has the biggest influence. By year-end, this year’s data takes over.

## Reading the Screen

The Simulator tab includes a few main sections.

* **The slider.** Lets you jump to any point in the year, or you can press Play to watch the year unfold month by month.
* **The 25-day breach rate.**  Shows how the model performed in the most recent month.
* **The year-to-date rate.** Displays breaches so far and shows today’s zone.
* **The true breach rate.** The app’s best estimate of how often the model actually fails. A good model should be close to 1%. Below, you’ll see the chance that the real rate is worse than 1%.
* **The zone forecast.** Uses a colored bar; it's the most important part of the page. It shows the odds of ending the year in green, yellow, or red.
* **The uncertainty curve.**  Shows each possible year-end breach count as a bar. Taller bars mean higher likelihood. If the bars are spread out, the outcome is still uncertain.
* **The clustering check.** Looks for breaches that happen close together, which is a warning sign. In a healthy model, misses should be spread out. If they cluster, this flag turns on.

## A Story in Three Moments

The sample data shows a story that’s worth exploring.

**Before the year begins.** History looks good. The app gives an 84% chance of green.

**Day 150.** There have been 4 breaches, so the official score is still green. However, the app now gives only a 23% chance of finishing green and a 76% chance of yellow. The clustering flag is also on. This is the early warning: the scoreboard says green, but the forecast warns that trouble is coming.

**Day 175.** Three more breaches occur. Green is now out of reach. The model ends the year with 7 breaches, landing in yellow.

A board that only watched the traffic light would have been surprised at year-end. With this app, they would have known months in advance.

## Things to Try

* Press Play to watch the zone bar change month by month.
* Change the breaches in one month from 0 to 3 and see how fast the forecast responds.
* Open “Prior source & weights” and select “Model’s claim (1%).” This means the app starts by fully trusting the model. Notice how this slows the alarm but doesn’t stop it.

## Questions to Ask Your Risk Team

* What is the chance we end this year in yellow or red?
* Are our breaches spread out, or clustered together?
* If we land in yellow, how much more capital would we need to hold?
* How much of our confidence comes from past years rather than this year?

## Limits to Keep in Mind

* **This is a teaching tool.** It uses sample data, not actual results from a real bank.
* **Clustering makes things look better than they are.** The forecast assumes breaches happen independently. If they cluster, the real range of outcomes is wider than the app suggests.
* **The rules have changed.** The current Basel rules use the same breach count but a different penalty scale. This app uses the original 1996 version because it’s easier to learn from.
* **A forecast is not a verdict.** “76% chance of yellow” means it’s likely, but not certain.
