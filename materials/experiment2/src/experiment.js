// Experimental trials code
const timeline = [];

// This code is very repetitive. It will be good to change this in the future:
    const run_experiment = () => {
        // Setting parameters
        // This seems to be the solution to the problem with jatos
        // clors and blocks needs to be global
        const urlvar = (jatos_run) ? jatos.urlQueryParameters: jsPsych.data.urlVariables();
        const norew = (urlvar.phase != undefined && ["Reversal", "Devaluation", "Omission", "Extinction"].includes(capitalize(urlvar.phase))) ?
            capitalize(urlvar.phase) :
            "Reversal";
        const blocks = (Number(urlvar.blocks) == 0) ? 0 : (!isNaN(Number(urlvar.blocks))) ? Number(urlvar.blocks) : 12; 
        const prac = true; 
        const gam = true;


        if (urlvar.phase == undefined) console.log("No phase parameter used. Default is Reversal.")
        else if (!["Reversal", "Devaluation", "Omission", "Extinction"].includes(capitalize(urlvar.phase))) console.log(`WARNING: an invalid phase parameter was used: ${urlvar.phase}. Phase has been set to Reversal.`);

        console.log(`Experiment Parameters: Phase: ${norew}. Blocks: ${blocks}. Practice: ${prac}`);
        var experimentType = ""; 
        console.log("Experimento:", experimentType);

        trialObj = create_trials(blocks, norew, prac);
        console.log(trialObj)
        const [colorHigh, colorLow] = (blocks != 0) ? trialObj["Reward"][1].colors : ["orange", "blue"];
        if (blocks != 0) console.log(`Color high is ${colorHigh}. Color low is ${colorLow}.`)  


    
        // Experiment 2a or 2b
        const experiment = {
            type: jsPsychHtmlButtonResponse,
            stimulus: `
                <div style="margin-bottom: 50px;">
                    <h3>Elige la versión que quieres utilizar:</h3>
                </div>
            `,
            choices: [
                '<div style="margin-left: 10px; margin-right: 10px; margin-top: 20px; margin-bottom: 20px;"> Experimento<br><br><b>2 A</b> </div>', 
                '<div style="margin-right: 10px; margin-left: 10px; margin-top: 20px; margin-bottom: 20px;"> Experimento<br><br><b>2 B</b> </div>'
            ],
            //on_start: document.body.classList.add("gray"),
            on_finish: (data) => {
                const responseMapping = {
                    0: 'A', // Experiment 2a: nominal reversal
                    1: 'B'  // Experiment 2b: without instructions, just feedback
                };
                
                experimentType = responseMapping[data.response],
                jsPsych.data.addProperties({
                    Experiment: responseMapping[data.response],
                });
                console.log("Experimento:", experimentType);
            }
        };

        // University
        const university = {
            type: jsPsychHtmlButtonResponse,
            stimulus: `
                <div style="margin-bottom: 50px;">
                    <h3>¿En qué universidad estás?</h3>
                </div>
            `,
            choices: ['Universidad de Málaga <br><br> UMA', 'Universidad de Granada <br><br> UGR'],
            on_finish: (data) => {
                jsPsych.data.addProperties({
                    University: data.response, // University = 0 (Málaga); = 1 (Granada)
                })
            }
        };

        /*const codigoParticipante = {
            type: jsPsychHtmlButtonResponse,
            stimulus: function() {
                const generatedCode = generateRandomCode(); 
                const codeElement = `<div style="margin-bottom: 40px">Este es tu código de participante. Por favor, apúntalo antes de continuar:</div>
                <div style="border: 2px solid #ccc; padding: 50px; margin-bottom: 100px; width: 300px; margin-left: auto; margin-right: auto;">
                <p style="font-size: 1.5rem; font-weight: bold; margin: 0;">${generatedCode}</p></div>`
                jsPsych.data.addProperties({Code: generatedCode})
                //<p style="display: block; margin-bottom: 50px">Pulsa ${(lab) ? `<b>continuar</b> para comenzar el experimento` : `<b>terminar</b> para salir del experimento`}</p>`;
                return codeElement;
            },
            choices: ['Continuar'],
            on_finish: (data) => {
                if (jatos_run) {
                    const results = jsPsych.data.get().filter([{ trial_type: "psychophysics" }, { trial_type: "survey-html-form" }]).json();
                    jatos.submitResultData(results);
                }
            }
        };      */  


        // Participant ID:
        console.log("ID: " + urlvar.ID);
        jsPsych.data.addProperties({
            ID: urlvar.ID || "not provided",
        });


        // Experimental trials

        // Points necessary for earning medals. If extinction, points cut-offs are not increased 
        const points_cut_off = [15500, 27200, 31900, 37300, 40000, 49300, 57000].map((p) => {
            return (norew == "Reversal" & p != 15500) ? p * 1.3 : p; 
        })
        
        let trialNum = 0, total_points = 0, BlockNum = 0;

        const trial = {
            type: jsPsychPsychophysics,
            stimuli: () => {
                const sF = (lab)? 40: jsPsych.data.get().last(1).values()[0].px2deg;
                const log = jsPsych.timelineVariable("trialLog");
                console.log(log);
                // Stimulus size is determined to an scaling factor that transform pixels to degrees of visual angle
                return draw_display(1.15 * sF, 0.2 * sF, 5.05 * sF, log, jsPsych.timelineVariable("colors"), jsPsych.timelineVariable("orientation"));
            },
            choices: ['g', 'c'],
            background_color: '#000000',
            canvas_width: () => {
                const sF = (lab)? 40: jsPsych.data.get().last(1).values()[0].px2deg;
                return sF * 15;
            },
            canvas_height: () => {
                const sF = (lab)? 40: jsPsych.data.get().last(1).values()[0].px2deg;
                return sF * 15;
            },
            data: () => {
                let color
                if (jsPsych.timelineVariable("condition") == "Absent") {
                    color = "none";
                } else {
                    if (jsPsych.timelineVariable("Phase") != "Reversal") {
                        color = (jsPsych.timelineVariable("condition") == "High") ? colorHigh : colorLow;
                    } else {
                        color = (jsPsych.timelineVariable("condition") == "Low") ? colorLow : colorHigh;
                    }
                }
                return {
                    tPos: jsPsych.timelineVariable("targetPos"),
                    sPos: jsPsych.timelineVariable("singPos"),
                    Phase: jsPsych.timelineVariable("Phase"),
                    condition: jsPsych.timelineVariable("condition"),
                    Block_num: (trialNum % 24 == 0) ? ++BlockNum : BlockNum, 
                    trial_num: ++trialNum,
                    counterbalance: counterbalance,
                    color: color,
                }
            },
            on_finish: (data) => {
                data.correct_response = (jsPsych.timelineVariable("orientation") == "vertical") ? "g" : "c";
                data.correct = (jsPsych.pluginAPI.compareKeys(data.response, data.correct_response)) ? 1 : 0;
                data.points = (data.correct) ?
                    compute_points(data.rt, data.condition, data.Phase) : -compute_points(data.rt, data.condition, data.Phase); //Con castigo: resta los mismos puntos que computa
                total_points = total_points + data.points; //(total_points + data.points <= 0) ? 0 : total_points + data.points;
                data.total_points = total_points;
            },
            trial_duration: 3700,
            response_start_time: 1700, 
        };


        const feedback = {
            type: jsPsychHtmlKeyboardResponse,
            stimulus: () => {
                const response = jsPsych.data.get().last(1).values()[0].key_press;
                if (response !== null) {
                    const acc = jsPsych.data.get().last(1).values()[0].correct;

                    if (jsPsych.timelineVariable("Phase") == "Practice" || jsPsych.timelineVariable("Phase") == "Omission" || jsPsych.timelineVariable("Phase") == "Reversal" && jsPsych.data.get().last(1).values()[0].Experiment == "A") {
                        return (acc) ? `<p style="color: yellow; font-size: 2rem;">Correcto</p>` :
                            `<p style="color: red; font-size: 2rem;">Error</p>`;
                    } // This part is running while NOMINAL REVERSAL phase

                    const bonus = (
                        (jsPsych.timelineVariable("condition") == "High" &&
                        (jsPsych.timelineVariable("Phase") != "Extinction" && jsPsych.timelineVariable("Phase") != "Devaluation" && jsPsych.timelineVariable("Phase") != "Reversal")) ||
                        (jsPsych.timelineVariable("condition") == "Low" && jsPsych.timelineVariable("Phase") == "Reversal" && jsPsych.data.get().last(1).values()[0].Experiment == "B")
                    ) ?
                        `<div style="background-color: ${(acc) ? `yellow` : `red`}; color: ${(acc) ? `black` : `black`}; font-size: 2rem; font-weight: ${(acc) ? `600` : `600`}; padding: ${(acc) ? `40px` : `40px`};">${(acc) ? //padding sirve para añadir espacios alrededor del texto 
                            `¡Puntos Extra!` :
                            `Pérdidas Extra`}</div><br>` :
                        '<div></div></br>';
                    const points = jsPsych.data.get().last(1).values()[0].points;
                    const gains = (acc) ?
                    `<p style="color: yellow; font-size: 2rem;">+${points} puntos</p>` :
                    `<p style="color: red; font-size: 2rem;">${points} puntos</p>`;

                    if ((jsPsych.timelineVariable("condition") != "High" && jsPsych.timelineVariable("Phase") == "Reward") || (jsPsych.timelineVariable("Phase") == "Extinction") || (jsPsych.timelineVariable("condition") != "Low" && jsPsych.timelineVariable("Phase") == "Reversal") && jsPsych.data.get().last(1).values()[0].Experiment == "B") {
                        const message = (acc) ?
                        `<p style="color: yellow; font-size: 2rem;">Correcto</p><br>` :
                            `<p style="color: red; font-size: 2rem;">Error</p><br>`;
                        return message + gains; // This part is running while Rewarded phase
                    } 
                    return bonus + gains;
                }
                return "<p style='font-size: 2rem;'>Demasiado lento. Intenta responder más rápido.</p>"
            },
            post_trial_gap: () => {
                const phase = jsPsych.data.get().last(1).values()[0].Phase;
                if ((phase == "Practice" && trialNum == 24) || (phase != "Reversal" && trialNum == (24 * blocks) * 2)) return 1000
            },
            trial_duration: 700,
            choices: ["NO_KEYS"],
        }

        const rest = {
            type: jsPsychHtmlKeyboardResponse,
            stimulus: () => {
                const rank = find_ranking(points_cut_off, total_points);
                const disp_medals = (rank >= 0) ? `
                            <p>Has desbloqueado la siguiente medalla:</p>
                            <img src="src/img/medals/${get_medal(rank)}" width="150" height="150">` : "";
                const next = (rank <= 4) ?
                    `<p>Te quedan ${formatting(next_points(points_cut_off, rank + 1, total_points).toString())} puntos para desbloquear la siguiente medalla.</p>` :
                    "";

                return `
                        <p>Has completado ${BlockNum} de ${blocks * 2} Bloques.</p>
                        <p>Llevas ${formatting(total_points.toString())} puntos acumulados.</p>
                        ${disp_medals + next}
                        <p>Pulsa la barra espaciadora cuando quieras continuar.</p>
                        `
            },
            choices: [' '],
            on_finish: () => {
                if (jatos_run) {
                    const results = jsPsych.data.get().filter([{ trial_type: "psychophysics" }, { trial_type: "survey-html-form" }]).json();
                    jatos.submitResultData(results);
                }
            }
        }

        const transition = {
            type: jsPsychHtmlKeyboardResponse,
            // Hay que escalar los puntos en función de la segunda fase
            stimulus: () => {
                if (norew == "Reversal") {
                    return `<p>Has terminado la primera mitad del experimento.</p>
                            <p>A partir de ahora van a cambiar las reglas que determinan los puntos que puedes ganar:</p>
                            <p>Si aparece el color ${colors_t(colorLow)}, ganarás 10 veces más puntos.</p>
                            <p>Si aparece el color ${colors_t(colorHigh)}, no ganarás puntos extra.</p>
                            <p>Además, no podrás ver en pantalla los puntos que ganes (solo aparecerá correcto o incorrecto), pero ¡seguirás ganando puntos! Podrás comprobar los puntos totales al final del bloque, como en la fase anterior.</p>
                            <p>Pulsa la barra espaciadora para continuar con el experimento.</p>`
                }
                if (norew == "Extinction") {
                    return `<p>Has terminado la primera fase del experimento.</p>
                            <p style = "width: 80%; margin: auto">En la segunda fase vas a continuar con la misma tarea (determinar la orientación de la línea dentro del rombo), pero a partir de ahora el círculo <b>${colors_t(colorHigh)}</b> ya no señala la disponibilidad de puntos extra. 
                            Ganarás o perderás los puntos que correspondan a tu rapidez en aciertos, 
                            como en el caso de los ensayos donde no hay colores o en los que aparece el círculo <b>${colors_t(colorLow)}</b>. En resumen, todo es igual salvo que el círculo <b>${colors_t(colorHigh)}</b> ya no significa que los puntos ganados o perdidos se multipliquen por diez.</p>
                            <p>Pulsa la barra espaciadora para continuar con el experimento.</p>`
                } if (norew == "Omission") {
                    return `<p>Has terminado la primera mitad del experimento.</p>
                            <p>Ahora vas a continuar con la tarea, <b style = "color:red"> pero ya no ganarás más puntos</b>.</p>
                            <p>Pulsa la barra espaciadora para continuar.</p>`
                }
                return `<p>Has terminado la primera mitad del experimento.</p>
                        <p>Ahora vas a continuar con la tarea, pero si aparece el círculo de ${colors_t(colorHigh)} la cantidad de puntos que ganarás será siempre 0.</p>
                        <p>Pulsa la barra espaciadora para continuar.</p>`;
            },
            choices: [" "],
            on_start: () => {document.body.classList.add("black")},
            on_finish: () => {
                if (jatos_run) {
                    const results = jsPsych.data.get().filter([{ trial_type: "psychophysics" }, { trial_type: "survey-html-form" }]).json();
                    jatos.submitResultData(results);
                }
            },
            conditional_function: () => {
                return (jsPsych.data.get().last(1).values()[0].Experiment != "B");
            },
        }

        const report = {
            type: jsPsychHtmlKeyboardResponse,
            stimulus: () => {
                const rank = find_ranking(points_cut_off, total_points, r = true);
                const medal = (rank >= 0) ? `<p>Has acumulado ${formatting(total_points.toString())} puntos y desbloqueado la siguiente medalla: </p>
                        <img src="src/img/medals/${get_medal((rank > 5) ? rank - 1 : rank)}" width="150" height="150">` :
                    `<p>Has acumulado ${formatting(total_points.toString())} puntos.</p>`;
                return `<p>Acabas de terminar esta tarea.</p>
                        ${medal + report_performance(rank)}
                        <p>Pulsa la barra espaciadora para seguir.</p>`
            },
            choices: [' '],
            on_finish: () => {
                if (jatos_run) {
                    const results = jsPsych.data.get().filter([{ trial_type: "psychophysics" }, { trial_type: "survey-html-form" }]).json();
                    jatos.submitResultData(results);
                }
                document.body.classList.remove("black");
                document.body.classList.add("gray");
                document.body.style.cursor = 'auto';
            },
            post_trial_gap: 1000,
        };

        // Awareness test
        
        const trial_awareness = {
            type: jsPsychPsychophysics,
            stimuli: () => {
                const sF = (lab)? 40: jsPsych.data.get().last(1).values()[0].px2deg;
                const log = jsPsych.timelineVariable("trialLog");
                console.log(log);
                // Stimulus size is determined to an scaling factor that transform pixels to degrees of visual angle
                return draw_display(1.15 * sF, 0.2 * sF, 5.05 * sF, log, jsPsych.timelineVariable("colors"), jsPsych.timelineVariable("orientation"));
            },
            choices: ['c', 'm'],
            background_color: '#000000',
            canvas_width: () => {
                const sF = (lab)? 40: jsPsych.data.get().last(1).values()[0].px2deg;
                return sF * 15;
            },
            canvas_height: () => {
                const sF = (lab)? 40: jsPsych.data.get().last(1).values()[0].px2deg;
                return sF * 15;
            },
            data: () => {
                let color
                if (jsPsych.timelineVariable("condition") == "Absent") {
                    color = "none";
                } else {
                    if (jsPsych.timelineVariable("Phase") != "Reversal") {
                        color = (jsPsych.timelineVariable("condition") == "High") ? colorHigh : colorLow;
                    } else {
                        color = (jsPsych.timelineVariable("condition") == "Low") ? colorLow : colorHigh;
                    }
                }
                return {
                    tPos: jsPsych.timelineVariable("targetPos"),
                    sPos: jsPsych.timelineVariable("singPos"),
                    Phase: jsPsych.timelineVariable("Phase"),
                    condition: jsPsych.timelineVariable("condition"),
                    Block_num: (trialNum % 24 == 0) ? ++BlockNum : BlockNum, 
                    trial_num: ++trialNum,
                    counterbalance: counterbalance,
                    color: color,
                }
            },
            on_finish: (data) => {
                data.correct_response = (jsPsych.timelineVariable("condition") == "Low") ? "c" : "m"; // !!
                data.correct = (jsPsych.pluginAPI.compareKeys(data.response, data.correct_response)) ? 1 : 0;
            },
            //trial_duration: 3700,
            response_start_time: 1700, 
        };

        // Conditional 

        const if_call = {
            timeline: [call_experimenter],
            conditional_function: () => {
                return lab;
            },
        }

        const if_nodeRest = {
            timeline: [rest],
            conditional_function: () => {
                return (trialNum % 24 == 0 && trialNum != 24 * blocks) && (trialNum % 24 == 0 && trialNum != (24 * blocks) * 2);
            },
        }

        const rest_aw = {
            timeline: [
                {
                    type: jsPsychHtmlKeyboardResponse,
                    stimulus: ["<p>Has completado el primer bloque.</p><br><p>Pulsa la barra espaciadora cuando quieras continuar.</p>"],
                    choices: [' '],
                    on_finish: () => {
                        if (jatos_run) {
                            const results = jsPsych.data.get().filter([{ trial_type: "psychophysics" }, { trial_type: "survey-html-form" }]).json();
                            jatos.submitResultData(results);
                        }
                    }
                },
            ],
            conditional_function: () => {
                return (trialNum % 24 == 0 && trialNum != 24 * blocks)  && (trialNum % 24 == 0 && trialNum != (24 * blocks) * 2) && (trialNum % 24 == 0 && trialNum != ((24 * blocks) * 2) + 24 * 2);
            },
        }

        const reward = {
            timeline: [trial, feedback, if_nodeRest],
            timeline_variables: (blocks != 0) ? trialObj.Reward : [],
            repetitions: 1,
            randomize_order: false,
        }

        const noreward = {
            timeline: [trial, feedback, if_nodeRest],
            timeline_variables: (blocks != 0) ? trialObj[norew] : [],
            repetitions: 1,
            randomize_order: false,
        }

        const awareness_test = {
            timeline: [trial_awareness, rest_aw],
            timeline_variables: (blocks != 0) ? trialObj.Awareness : [],
            repetitions: 1,
            randomize_order: false,
        }

        const if_reward = {
            timeline: [reward],
            conditional_function: () => {
                return blocks != 0;
            }
        }

        const if_noreward = {
            timeline: [noreward],
            conditional_function: () => {
                return blocks != 0;
            }
        }

        const practice = {
            timeline: [trial, feedback],
            timeline_variables: trialObj["Practice"],
            repetitions: 1,
            randomize_order: false,
        }

        const if_practice = {
            timeline: [practice],
            conditional_function: () => {
                return prac;
            }
        }

        const procedure_cal = {
            timeline: [welcome, instructions_cal, full_on, resize],
            repetitions: 1,
            randomize_order: false,
        }

        const procedure_prac = {
            timeline: [instructions_prac, pre_prac, if_practice],
            repetitions: 1,
            randomize_order: false,
            post_trial_gap: () => { if (prac == false) return 1000 },
            on_finish: () => {
                console.log(trialNum);
                if (trialNum == 24 || prac == false) {
                    document.body.classList.remove("black");
                    document.body.classList.add("gray");
                    document.body.style.cursor = 'auto';
                    trialNum = 0;
                    BlockNum = 0;
                }
            }
        }

        const if_expA = {
            timeline: [instructions_exp, pre_exp, if_reward, if_call, transition, if_noreward],
            conditional_function: () => {
                return experimentType == "A"
            },
        }

        const if_expB = {
            timeline: [instructions_exp, pre_exp, if_reward, rest, if_noreward],
            conditional_function: () => {
                return experimentType == "B"
            },
        }

        const procedure_aw = {
            timeline: [start_awareness, instructions_aw, pre_aw, awareness_test],
            repetitions: 1,
            randomize_order: false,
            //post_trial_gap: 1000,
        }

        const download = {
            type: jsPsychHtmlButtonResponse,
            stimulus: `<p>¿Quieres descargar los datos?</p>`,
            choices: ["Sí", "No"],
            on_finish: (data) => {
                if (data.response == 0) {
                    jsPsych.data.get().filter([{ trial_type: "psychophysics" }, { trial_type: "survey-html-form" }]).localSave("csv", "example.csv");
                }
            }
        }

        const if_download = {
            timeline: [download],
            conditional_function: () => {
                return !jatos_run;
            }
        }


        const out_id = {
            type: jsPsychHtmlKeyboardResponse,
            stimulus: () => {
                //const id = codigoParticipante.response;
                //<p>Antes de terminar del experimento, apunta tu código de participante: ${id}</p>
                const html = (jsPsych.data.get().last(1).values()[0].University == 0) ? "<div>Has completado esta parte del experimento.</div> <br> <div> Permanece a la espera en tu sitio sin hacer nada. En breve, el experimentador te dará nuevas instrucciones.</div>" : "<p>Por favor, avisa al experimentador o experimentadora</p>"; 
                return html;
            },
            choices: ['s', 'q'],
            on_finish: (data) => {
                if (window.jatos && data.response == 'q') {
                    const results = jsPsych.data.get().filter([{ trial_type: "psychophysics" }, { trial_type: "survey-html-form" }]).json();
                    jatos.submitResultData(results)
                        .then(jatos.endStudyAjax)
                        .then(() => { window.location.href = 'https://encuestas.ugr.es/index.php/642174?lang=es&ID=' + jsPsych.data.get().last(1).values()[0].ID}) //Enlace a Qs, añade ID a la URL
                        .catch(() => console.log("Something went wrong"));    
                }
            }
        }
        
        const slide = {
            type: jsPsychHtmlSliderResponse,
            stimulus: () => { //YA NO TIENE MUCHO SENTIDO
                return `<div style="width:auto; margin-bottom: 50px;">
                <p>¿Qué porcentaje de puntos crees que has ganado con cada color?</p>
                <div style="width:240px; float: left;">
                    <canvas id="myCanvas1" width="150" height="150" style = "border-radius: 3%; background-color: #fff; margin: 10px 0;"></canvas>
                </div>
                <div style="width:240px; float: right;">
                    <canvas id="myCanvas2" width="150" height="150" style = "border-radius: 3%; background-color: #fff; margin: 10px 0;"></canvas>
                </div>
                </div>`
            },
            require_movement: true,
            labels: [`<span id="high-placeholder">50</span> % con el color ${colors_t(colorHigh)}`,
                `<span id="low-placeholder">50</span> % con el color ${colors_t(colorLow)}`],
            prompt: "<p>Pulsa continuar cuando hayas acabado</p>",
            button_label: "Continuar",
                on_load: () => {
                document.addEventListener("click", slider_c);
                const slider = document.getElementsByClassName("jspsych-slider");
                slider[0].addEventListener("input", slider_move);

            },
            on_finish: (data) => {
                document.removeEventListener("click", slider_c);
                jsPsych.data.addProperties({
                    contingency_rating: data.response,
                })
            },
            //post_trial_gap: 500,
        };

        const slider_proc = {
            timeline: [slider_instr, slide],
        }

        timeline.push(check, experiment, university, welcome, procedure_prac, if_call, if_expA, if_expB, report, procedure_aw, full_off, questions, out_id, if_download);
        jsPsych.run(timeline);
    }

    if (jatos_run) {
        jatos.onLoad(() => run_experiment());
    } else {
        run_experiment();
    }
