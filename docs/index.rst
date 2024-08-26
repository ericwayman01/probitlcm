=====================================
Welcome to probitlcm's documentation!
=====================================

This software implements the model described in the `manuscript <http://arxiv.org/abs/2408.13143>`_

  Eric Alan Wayman, Steven Andrew Culpepper, Jeff Douglas, and Jesse Bowers. "A Restricted Latent Class Model with Polytomous Ordinal Correlated Attributes and Respondent-Level Covariates." arXiv preprint arXiv:2408.13143, 2024.

If you make use of this software, please cite the above.

Introduction
============

The simulations and data analyses should be run from "run directories" which contain the necessary input files for the run. There are two types of run directories, one for simulations and one for data analyses. They are referred to as ``run_directory`` and ``run_dir_data_analysis``, respectively.

Simulation
==========

Creating directories and files
------------------------------

First, the necessary directories and folders must be created, as follows.

Contents of ``run_directory`` directory:

- file: ``config_simulation.toml``
- directory: ``<sim_info_dir_name>``

::
   
   # My config file
   # these should be absolute paths
   # current run_directory used is
   #     laptop: /home/username/run_directory
   #     cluster: /home/username/run_directory
   
   laptop_process_dir = "/home/username/process_dir"
   cluster_process_dir = "/scratch/users/username/process_dir"
   number_of_replics = 100

Each simulation has its own directory, ``<sim_info_dir_name>``.

Contents of ``<sim_info_dir_name>``:

- ``json_files`` directory

Contents of ``json_files`` directory:

- ``01_fixed_vals.json``
- ``01_list_higherlevel_datagen_vals.json``
- ``01_list_hyperparams_possible_vals.json``
- ``01_list_N_vals.json``
- ``02_current_tuning_hyperparam_vals.json``
- ``02_list_decided_hyperparam_vals.json``

Example ``01_fixed_vals.json``:

::
   
   {
       "T": 1,
       "covariates": "age_assignedsex",
       "chain_length_after_burnin": 10000,
       "burnin": 6000,
       "q": 5,
       "order": 2,
       "omega_0": 0.5,
       "omega_1": 0.5,
       "sigma_beta_sq": 2.0
   }
   
Example ``01_list_higherlevel_datagen_vals.json``:

::
   
   {
       "L_k_s": [[2, 2], [3, 3],
                 [2, 2, 2], [3, 3, 3],
                 [2, 2, 2, 2]],
       "rho": [0, 0.25, 0.5]
   }
   

Example ``01_list_hyperparams_possible_vals.json``:

::
   
   {
       "sigma_kappa_sq": [0.00048828125, 0.0009765625, 0.001953125, 0.00390625, 0.0078125, 0.015625, 0.03125, 0.0625, 0.125]
   }
   

Example ``01_list_N_vals.json``:

::
   
   [500, 1500, 3000]
   
Example ``02_current_tuning_hyperparam_vals.json``:

::
   
   {
       "sigma_kappa_sq": 0.125
   }
   
Example ``02_list_decided_hyperparam_vals.json``:

::
   
   {
       "sigma_kappa_sq": [
           [
               500,
               1500,
               3000
           ],
           [
               0.00390625,
               0.001953125,
               0.0009765625
           ]
       ]
   }

Creating data generating parameters and other necessary files
-------------------------------------------------------------
   
Once the above directories and files have been created, the ``probitlcm_extra`` package (available at https://github.com/ericwayman01/probitlcm_extra) will be used to create the data generating parameters and other necessary files. Install that package, and then continue as follows:

From ``run_directory``, run, for example:

::
   
   python3.X -m probitlcm_extra.build_situations_list --sim_info_dir sim_for_manuscript

where ``sim_for_manuscript`` is an example of ``<sim_info_dir_name>``. ``02_list_situations.json`` will be placed in ``<run_directory>/<sim_info_dir_name>/json_files`` as a result.

Now, from ``run_directory``, for a particular ``sim_info_dir_name``, for ``situation_num`` ranging from 1 to the total number of situations (a situation is a particular combination of user-selected non-random values for a scenario to be simulated other than sample size; there are 15 situations in the manuscript).

::

   python3.X -m probitlcm_extra.generate_lambda --sim_info_dir_name sim_for_manuscript --situation_num 1

Finally run the command, for example

::

   python3.X -m probitlcm_extra.create_scenario_params --sim_info_dir sim_for_manuscript --situation_num 1

for each situation to create the necessary data-generating parameters.
   
Performing hyperparameter tuning
--------------------------------

Once one has prepared the files for the various scenarios (situations combined with sample sizes), before running simulations, one can perform hyperparameter tuning on a particular scenario. The list of hyperparameters and their possible values should be placed in ``<run_directory>/<sim_info_dir_name>/json_files/01_list_hyperparams_possible_vals.json``. For example:

::

   {
       "sigma_kappa_sq": [0.50, 0.75, 1.0, 1.5, 2.0]
   }
   
One runs

::

   python3.X -m probitlcm.hyperparam_tuning --environ laptop --dir_name sim_for_manuscript --type_of_run simulation --scenario_or_setup_number 1

The output will appear in the processing directory, in ``<process_dir>/<sim_info_dir_name>/hyperparam_tuning/scenario_<scenarionum>_<hyperparamname>_results.txt``. Each line of the file contains the following content:

::

   <hyperparam value>,<percentage of acceptances/rejections>

Currently, for the tuning, the software uses a "warmup" period of 1000 iterations followed a period of 2000 iterations for which the data forming the acceptance percentage is collected.

Based on this, the user should manually add the decided-upon value(s) of the hyperparameter(s) to the ``02_list_decided_hyperparam_vals.json`` file.

::

   {
       "sigma_kappa_sq": [
           [
               500,
               1500,
               3000
           ],
           [
               0.00390625,
               0.001953125,
               0.0009765625
           ]
       ]
   }

Running a simulation
--------------------

The aforementioned files and directories must have been created in order for one to run a simulation.

A scenario is run in two pieces: the ``scenario_launch`` function launches all replications for a given scenario. In a laptop environment, control is returned to the user following the completion of all replications which are run one-by-one. In a cluster environment, all replications are launched simultaneously and control is then returned to the user.

Once the user thinks the replications are completed, the user runs the ``scenario_replics_postprocess`` function to perform post-processing.

Some example simulation run commands:

::
   
   python3.X -m probitlcm.scenario_launch --environ laptop --sim_info_dir_name sim_for_manuscript --scenarionumber 1
   python3.X -m probitlcm.scenario_replics_postprocess --environ laptop --sim_info_dir_name sim_for_manuscript --scenarionumber 1

Once all scenarios for a simulation have finished running, a report on all scenarios can be built using the ``build_report_all_scenarios`` module. The command is

::

   python3.X -m probitlcm.build_report_all_scenarios --environ laptop --sim_info_dir_name sim_for_manuscript --total_num_of_scenarios 45

Data analysis
=============

Similarly to a simulation, the ``run_dir_data_analysis`` must exist, and have a configuration file as well as a directory for the data set for which one wishes to perform an analysis.

Contents of run directory
-------------------------

Contents of ``run_dir_data_analysis`` directory:

- file: ``config_data_analysis.toml``
- directory: ``<dataset_name>``

Example ``config_data_analysis.toml``:
  
::
   
   # My config file
   # these should be absolute paths
   # current run_directory used is
   #     laptop: /home/username/run_dir_data_analysis
   #     cluster: /home/username/run_dir_data_analysis
   
   laptop_process_dir = "/home/username/process_dir_data_analysis"
   cluster_process_dir = "/scratch/users/username/process_dir_data_analysis"
   number_of_chains = 1
   
Contents of ``<dataset_name>`` directory:

- ``01_fixed_vals.json``
- ``01_list_hyperparams_possible_vals.json``
- ``02_decided_hyperparam_vals.json``
- ``covariates.csv``
- ``M_j_s.json``
- ``responses.csv``
- ``setup_0001.json``, ..., ``setup_000X.json``

Example ``01_fixed_vals.json``:

::
   
   {
       "T": 1,
       "covariates": "age_assignedsex",
       "chain_length_after_burnin": 60000,
       "burnin": 60000,
       "order": 2,
       "omega_0": 0.5,
       "omega_1": 0.5,
       "sigma_beta_sq": 2.0
   }

Example ``01_list_hyperparams_possible_vals.json``:

::

   {
    "sigma_kappa_sq": [0.00048828125, 0.0009765625, 0.001953125, 0.00390625, 0.0078125, 0.015625, 0.03125, 0.0625, 0.125]
   }

Example ``02_decided_hyperparam_vals.json``:

::

   {
       "sigma_kappa_sq": 0.001953125
   }

Example ``M_j_s.json``:

::

   [3, 3, 3, 5, 5, 3, 3, 4, 5, 5, 5, 5, 5, 3, 5, 5, 3]

Example ``setup_000x.json``:

::

   {
       "K": 2,
       "L_k_s": [
           2,
           2
       ]
   }
   
Explanation of contents
-----------------------
   
A dataset may have multiple setups (like scenarios for the simulations).

The ``.csv`` files must be comma-separated and should have a header. It is convenient to save such files from ``R`` using, for example, the function call ``write.csv(df, "nan_test.csv", row.names=FALSE)`` (an equivalent function should exist from Python but has not been tested). The option ``csv_opts::strict`` in Armadillo specifies that "missing" values will be converted to ``arma::datum::nan`` type and that seems to apply to the ``NA`` values (really, the strings ``"NA"``) placed in those cells by R. Presumably any non-numeric type would work, but such situations have not been tested.

The data in ``responses.csv`` must be zero-based, i.e. a column ``j`` must run from ``0`` through ``M_j - 1``.

``M_j_s.json`` should consist of a list whose length is the number of items in ``responses.csv`` and whose elements indicate the number of possible levels for each item. For example:

::
   
   [3, 3, 3, 5, 5, 3, 3, 4, 5, 5, 5, 5, 5, 3, 5, 5, 3]

Performing hyperparameter tuning
--------------------------------

One runs

::

   python3.X -m probitlcm.hyperparam_tuning --environ laptop --dir_name sim_for_manuscript --type_of_run data_analysis --scenario_or_setup_number 1

The content of the file ``02_decided_hyperparam_vals.json`` should be filled in based on the results of hyperparameter tuning.
   
Running a data analysis
-----------------------

An example data analysis run command:

::
   
   python3.X -m probitlcm.data_analysis_run --environ laptop --dataset rush_cross_sec --setup_num 1

Performing posterior predictive checks
--------------------------------------

Example command:

::
   
   python3.X -m probitlcm.posterior_predictive_check --environ cluster --dataset rush_cross_sec --setup_num 1 --chainnum 1 --N_files 6

The value for ``N_files`` must match the number of files containing draws from the posterior predictive distribution, which were created during the run of the particular ``setup_name``.

The run of ``posterior_predictive_check`` produces files ``bartholomew_histogram.png`` and ``bartholomew_p_value.txt`` in the particular chain's subdirectory.

Seed numbers
============

::
   
   seed_number = n_replics * (scenario_number - 1) + (replic_num - 1)

License
=======

``probitlcm`` is proved under the GNU General Public License v3, a copy of which can be found in the ``LICENSE`` file. By using, distributing, or contributing to this project, you agree to the terms and conditions of this license.

API Documentation
=================

.. toctree::
   :maxdepth: 2

   _core
   py_files

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
