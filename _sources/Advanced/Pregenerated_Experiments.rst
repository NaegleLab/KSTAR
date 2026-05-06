
Pregenerated Random Activities
===============================

In KSTARv1.0, we introduced the use of random experiments that were pre-generated and reused for future activity predictions, rather than generated at the time of activity prediction. This significantly reduced computational time, but does mean that the random experiments likely do not exactly match the true experiment as they did in the original implementation of KSTAR. 

By default, KSTARv1.0 will use pregenerated random experiments with a similar size and bias distribution as the real experiment (greater than 150 sites, within 25% of the pregenerated number of sites). If there are no appropriate pregenerated experiments, the random experiments will be generated from scratch.

Adjusting Pregenerated Experiments Use 
--------------------------------------

While we did not find significant impact on results when using pregenerated experiments, users may wish to default to the original, published version of the algorithm where random experiments are generated during runtime to match the real experiment exactly. Should you wish to adjust the default behavior of KSTAR, you can tweak the below parameters that can be passed to :func:`run_kstar_analysis()<kstar.calculate.run_kstar_analysis>`:

- `use_pregen_data`: set to False to not use pregenerated random experiments (default = True)
- `max_diff_from_pregenerated`: maximum fractional difference between the number of sites in the real and pregenerated experiments (default = 0.25)
- `min_dataset_size_for_pregenerated`: minimum number of sites required to use pregenerated random experiments
-`require_pregenerated`: if True, will ignore data columns that do not have an appropriate pregenerated experiment. This will help improve speed, but may cause dataset loss

Generating Additional Pregenerated Experiments
----------------------------------------------

In addition, you may wish to generate additional pre-generated experiments that can be used in the future for cases where there is not an appropriate pre-generated experiment already existing. This can be done two ways:

1. Automatically save newly generated random activities during KSTAR activity calculation based properties of the real samples
2. Manually generate a set of pregenerated random activities to be used in future KSTAR activity calculations with the pregenerate module

Automatically during KSTAR activity Calculation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The easiest way to create new experiments is to let KSTAR generate and save these automatically when there are no appropriate random experiments available. This can be done by tweaking the parameters provided in the :func:`run_kstar_analysis()<kstar.calculate.run_kstar_analysis>`, including indicating where you would like these custom activities stored:
- `custom_pregenerated_path`: where additional pre-generated experiments are stored and will be saved to if `save_new_random_activities=True`
- `save_new_random_activities`: save newly generated random experiments in the `custom_pregenerated_path`
- `default_pregen_only`: if True, will ignore any experiments in the `custom_pregenerated_path`

Manually with the pregenerate module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Some users may want a specific experiment set generated ahead of time. This can be done using the `kstar.pregenerate` module. Let's say you want to create a single random phosphotyrosine experiment with 500 sites and a study bias distribution of 5% low, 15% medium, and 80% high. You can do this with the following code::

    from kstar.random_experiments import pregenerate

    random_activities = pregenerate.generate_random_activities_from_scratch(
        exp_size = 500,
        compendia_sizes = (0.05, 0.15, 0.80),
        phospho_type = 'Y'
    )






