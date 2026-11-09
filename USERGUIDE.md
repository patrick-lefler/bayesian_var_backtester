# Bayesian VaR Backtester: A Plain-English Guide

## The Big Idea

Banks use a risk model called Value at Risk, or VaR, to answer one question: "How much could we lose on a bad day?" A 99% VaR is a promise. It says that on 99 out of 100 trading days, losses will stay below a set number.

On the other 1 day in 100, losses go past that number. That's called a **breach**. Breaches are expected. A good model should have about 2 or 3 in a year of 250 trading days.

Regulators keep score with a traffic light:

| Zone | Breaches in 250 days | What it means |
|---|---|---|
| Green | 0 to 4 | The model looks fine. |
| Yellow | 5 to 9 | Questions get asked, and the bank may have to hold more capital. |
| Red | 10 or more | The model is presumed broken, and the capital penalty is largest. |

The problem is that the score only counts at the end of the year. Halfway through, a count of 4 breaches looks "green," but that tells you little about where the year will finish. This app answers the question a board actually cares about: **"Where are we likely to end up, and how sure are we?"**

## How the App Thinks

The app works like a sensible person. It starts with what it already knows, then updates that view as new evidence arrives. This method is called Bayes' Theorem.

1. **Start with history.** The model had 3 breaches two years ago and 2 last year. Last year counts a little more because it's more recent.
2. **Add new evidence.** Every 25 trading days, about one month, the app adds the latest breaches to the record.
3. **Forecast the finish.** It works out the chances of ending the year green, yellow, or red.

Early in the year, history does most of the work. By the end, this year's data does.

## Reading the Screen

The Simulator tab has a few main parts.

- **The slider.** Move it to any point in the year, or press Play to watch the year unfold month by month.
- **25-day breach rate.** How the model did in the most recent month.
- **Year-to-date rate.** Breaches so far, and today's zone.
- **True breach rate.** The app's best guess at how often the model really fails. A sound model should be close to 1%. Underneath, it shows the chance that the real rate is worse than 1%.
- **Zone forecast.** The colored bar is the most important item on the page. It shows the odds of ending the year green, yellow, or red.
- **Uncertainty curve.** Each bar is one possible year-end breach count. Taller bars are more likely. A wide spread means the outcome is still open.
- **Clustering check.** Breaches that arrive in bunches are a warning sign, because a healthy model's misses should be scattered. If they bunch up, this flag lights up.

## A Story in Three Moments

The sample data tells a story worth walking through.

**Before the year starts.** History looks good. The app gives an 84% chance of green.

**Day 150.** There have been 4 breaches, so the official score is still green. But the app now gives only a 23% chance of finishing green and a 76% chance of yellow. The clustering flag is also on. **This is the early warning.** The scoreboard says green, but the forecast says trouble is coming.

**Day 175.** Three more breaches arrive. Green is now out of reach. The model ends the year with 7 breaches, in yellow.

A board that saw only the traffic light would have been surprised at year end. A board using this app would have known months earlier.

## Things to Try

- Press Play and watch the zone bar shift month by month.
- Change one month's breaches from 0 to 3 and see how quickly the forecast reacts.
- Open "Prior source & weights" and choose "Model's claim (1%)." This starts by trusting the model completely. Notice how it slows the alarm without stopping it.

## Questions to Ask Your Risk Team

- What is the chance we end this year in yellow or red?
- Are our breaches spread out, or clustered together?
- If we land in yellow, how much more capital would we need to hold?
- How much of our confidence comes from past years rather than this year?

## Limits to Keep in Mind

- **It's a teaching tool.** It uses sample data, not a real bank's results.
- **Clustering makes things worse than they look.** The forecast assumes breaches are independent. If they cluster, the real range of outcomes is wider than the app shows.
- **The rules have changed.** Current Basel rules use the same breach count but a different penalty scale. The app follows the original 1996 version because it's simpler to learn from.
- **A forecast isn't a verdict.** "76% chance of yellow" means likely, not certain.
