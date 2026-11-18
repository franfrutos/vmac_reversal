# Experiments 2a and 2b: nominal reversal and reversal without instructions

Experimental task of Experiments 2a and 2b of **"Explicit knowledge modulates Pavlovian attentional biases"**. It is the task of Experiment 1 (visual search, *additional singleton* paradigm, with an acquisition stage and a reversal stage), written in [jsPsych](https://www.jspsych.org/7.2/) 7.2.3 with the [jspsych-psychophysics](https://jspsychophysics.hes.kyushu-u.ac.jp/) plugin, with the changes of Experiments 2a and 2b and an awareness test at the end. Participants completed it in the laboratory at the University of Granada and the University of Málaga, and the task was served by the JATOS server of the University of Granada (May to October 2024).

At the start of each session, the experimenter selects the version of the task:

- **Experiment 2a (`A`)**: before the reversal stage, participants are informed that the color-reward contingencies are reversed. During the reversal stage, feedback only indicates whether the response was correct, and the points accumulated are shown at the end of each block.
- **Experiment 2b (`B`)**: participants are not informed about the reversal and keep receiving trial-by-trial feedback with points, which now follows the reversed contingencies.

## Quick start

- **Locally:** serve the folder (for example, `python -m http.server` in this folder) and open `http://localhost:8000/index.html?ID=123456`. At the end, the task offers to download the data as a CSV.
- **Online:** create a study in [JATOS](https://www.jatos.org/), copy the content of this folder into the assets folder of the study, and add a component whose HTML file is `index.html`. The data are sent to JATOS after each block and at the end of the session. The last screen can then redirect participants to a questionnaire on the survey server of the University of Granada, with their participant code in the link (`out_id` in `src/experiment.js`).

### URL parameters

| Parameter | Values | Default |
|---|---|---|
| `phase` | stage after acquisition: `Reversal`, `Devaluation`, `Omission` or `Extinction` | `Reversal` (the stage used in Experiments 2a and 2b) |
| `blocks` | number of blocks of each stage | 12 |
| `ID` | participant code | `not provided` |

## Settings

The main settings are in `src/init.js`:

| Setting | What it does |
|---|---|
| `lab` | `true` (laboratory): stimuli are drawn at 40 pixels per degree, and a screen asks participants to wait for the experimenter after the practice and, in Experiment 2a, before the reversal stage |
| `counterbalance` | pair of colors of the high- and low-value singletons, drawn at random for each participant |

The number of trials, the colors and the points are defined in `src/helper.js`. All texts shown to participants (in Spanish) are in `src/text.js`, and the instructions about the reversal of Experiment 2a are in the `transition` screen of `src/experiment.js`.

## Data

One row per search trial and awareness-test trial (`trial_type` = `psychophysics`) and one row for the post-experiment questionnaire (`survey-html-form`). The main columns are:

| Column | Content |
|---|---|
| `Experiment` | version of the task: `A` (Experiment 2a) or `B` (Experiment 2b) |
| `University` | `0`: University of Málaga; `1`: University of Granada |
| `Phase` | `Practice`, `Reward` (acquisition), `Reversal` or `Awareness` |
| `Block_num`, `trial_num` | block and trial |
| `Singleton`, `color` | singleton distractor (`High`, `Low` or `Absent`) and its color (`none` if absent) |
| `tPos`, `sPos` | position of the target and of the singleton (`-1`: absent) |
| `rt`, `key_press`, `correct_response`, `correct` | response time (ms), key pressed, correct key and accuracy |
| `points`, `total_points` | points of the trial and accumulated points |
| `distraction_rating`, `opinion_text` | answers to the post-experiment questions |

Every row also carries the participant code (`ID`), the random seed (`rng_seed`) and `counterbalance`.

## Folder structure

```
index.html          loads jsPsych, the plugins and the task
src/init.js         settings and initialization of jsPsych
src/helper.js       stimuli, trial lists, colors and points
src/text.js         texts shown to participants and text-only screens
src/experiment.js   the session: version, instructions, practice, acquisition, reversal, awareness test and questions
src/stylesheet.css  styles
src/img/            images of the awareness test (bonus.png, noBonus.png), medals shown between blocks, and the card of the
                    virtual chin-rest (dni.jpg), which this version does not use
jspsych/            jsPsych 7.2.3, the plugins loaded by index.html and jspsych-psychophysics
```
