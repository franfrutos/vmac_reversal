# Experiment 1: instructed reversal

Experimental task of Experiment 1 of **"Explicit knowledge modulates Pavlovian attentional biases"**. It is a visual search task (*additional singleton* paradigm) adapted from Le Pelley et al. (2015) and Watson et al. (2019), with an acquisition stage and a reversal stage, written in [jsPsych](https://www.jspsych.org/7.2/) 7.2.3 with the [jspsych-psychophysics](https://jspsychophysics.hes.kyushu-u.ac.jp/) plugin. Participants completed it in the laboratory, and the task was served by the JATOS server of the University of Granada.

## Quick start

The task runs only inside [JATOS](https://www.jatos.org/), because it starts when JATOS loads the page (`jatos.onLoad` in `src/experiment.js`). To run it, create a study in JATOS (a local installation of JATOS also works), copy the content of this folder into the assets folder of the study, and add a component whose HTML file is `index.html`. The data are sent to JATOS after each block and at the end of the session.

### URL parameters

| Parameter | Values | Default |
|---|---|---|
| `phase` | stage after acquisition: `Reversal`, `Devaluation` or `Omission` | `Reversal` (the stage used in Experiment 1) |
| `blocks` | number of blocks of each stage | 12 |
| `ID` | participant code | `not provided` |

## Settings

The main settings are in `src/init.js`:

| Setting | What it does |
|---|---|
| `lab` | `true` (laboratory): stimuli are drawn at 40 pixels per degree, and a screen asks participants to call the experimenter after the practice and before the reversal stage |
| `counterbalance` | pair of colors of the high- and low-value singletons, drawn at random for each participant |

The number of trials, the colors and the points are defined in `src/helper.js`. All texts shown to participants (in Spanish) are in `src/text.js`, and the instructions about the reversal are in the `transition` screen of `src/experiment.js`.

## Data

One row per search trial (`trial_type` = `psychophysics`) and one row for the post-experiment questionnaire (`survey-html-form`). The main columns are:

| Column | Content |
|---|---|
| `Phase` | `Practice`, `Reward` (acquisition) or `Reversal` |
| `Block_num`, `trial_num` | block and trial |
| `Singleton`, `color` | singleton distractor (`High`, `Low` or `Absent`) and its color |
| `tPos`, `sPos` | position of the target and of the singleton (`-1`: absent) |
| `rt`, `key_press`, `correct_response`, `correct` | response time (ms), key pressed, correct key and accuracy |
| `points`, `total_points` | points of the trial and accumulated points |
| `contingency_rating`, `distraction_rating`, `opinion_text` | answers to the post-experiment questions |

Every row also carries the session information: participant code (`ID`), random seed (`rng_seed`), `counterbalance`, and the calibration of the virtual chin-rest (`px2deg`, `viewing_distance`).

## Folder structure

```
index.html          loads jsPsych, the plugins and the task
src/init.js         settings and initialization of jsPsych
src/helper.js       stimuli, trial lists, colors and points
src/text.js         texts shown to participants and text-only screens
src/experiment.js   the session: calibration, instructions, practice, acquisition, reversal and questions
src/stylesheet.css  styles
src/img/            card for the calibration (dni.jpg), medals shown between blocks, and images of an unused memory task
jspsych/            jsPsych 7.2.3, the plugins loaded by index.html and jspsych-psychophysics
```

## References

- Le Pelley, M. E., Pearson, D., Griffiths, O., & Beesley, T. (2015). When goals conflict with values: Counterproductive attentional and oculomotor capture by reward-related stimuli. *Journal of Experimental Psychology: General*, *144*(1), 158–171. https://doi.org/10.1037/xge0000037
- Watson, P., Pearson, D., Most, S. B., Theeuwes, J., Wiers, R. W., & Le Pelley, M. E. (2019). Attentional capture by Pavlovian reward-signalling distractors in visual search persists when rewards are removed. *PLOS ONE*, *14*(12), e0226284. https://doi.org/10.1371/journal.pone.0226284
