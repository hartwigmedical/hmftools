package com.hartwig.hmftools.cup.runners;

import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;

import java.io.File;
import java.io.IOException;

import com.google.common.annotations.VisibleForTesting;

import org.apache.commons.io.FileUtils;
import org.apache.logging.log4j.Level;

public class PythonEnvironment
{
    public static final String DEFAULT_PYTHON_VERSION = "3.9.4";

    private static final File PYENV_DIR = new File(System.getProperty("user.home") + "/.pyenv");
    private static final File PYENV_PATH = new File(PYENV_DIR + "/libexec/pyenv");

    public final String mPythonVersion;
    public final File mVirtualEnvName;
    public final boolean mQuietInitialize;

    public PythonEnvironment(String pythonVersion, String virtualEnvName, boolean quietInitialize)
    {
        mVirtualEnvName = new File(virtualEnvName);
        mPythonVersion = pythonVersion;
        mQuietInitialize = quietInitialize;
    }

    public File virtualEnvPath() { return new File(String.format("%s/versions/%s/envs/%s", PYENV_DIR, mPythonVersion, mVirtualEnvName)); }

    public File pythonPath() { return new File(virtualEnvPath() + "/bin/python"); }

    @VisibleForTesting
    void installPyenv()
    {
        if(!PYENV_DIR.exists())
        {
            String command;

            command = "curl https://pyenv.run | bash"; // This also installs pyenv-virtualenv to $HOME/.pyenv/plugins/pyenv-virtualenv
            CUP_LOGGER.info("Installing pyenv with command: " + command);
            new BashCommand(command).logLevel(Level.DEBUG).run();
        }
        else if(!mQuietInitialize)
        {
            CUP_LOGGER.warn("Skipping installing pyenv as it exists at: " + PYENV_DIR);
        }

        if(!PYENV_DIR.exists())
        {
            CUP_LOGGER.error("Failed to install pyenv to " + PYENV_DIR);
            System.exit(1);
        }
    }

    @VisibleForTesting
    void addPyenvPathsToRcFile()
    {
        String shellPath = new BashCommand("echo $SHELL").logLevel(null).run().getStdout().get(0);
        String shellName = new File(shellPath).getName();
        File RcFilePath = new File(String.format("%s/.%src", System.getProperty("user.home"), shellName));

        if(!RcFilePath.isFile())
        {
            CUP_LOGGER.error("Fail to find rc file({}) to add pyenv paths. Using pyenv in interactive shell will not work", RcFilePath);
            System.exit(1);
        }

        try
        {
            String lines = String.join("\n",
                "# Pyenv paths",
                "export PYENV_ROOT="+ PYENV_DIR,
                "[[ -d $PYENV_ROOT/bin ]] && export PATH=\"$PYENV_ROOT/bin:$PATH\"",
                "eval \"$(pyenv init -)\"",
                "eval \"$(pyenv virtualenv-init -)\""
            );

            boolean linesExist = FileUtils.readFileToString(RcFilePath, "UTF-8").contains(lines);
            if(!linesExist)
            {
                String command;

                CUP_LOGGER.info("Adding pyenv paths to rc file: " + RcFilePath);
                command = String.format("echo -e '\n%s' >> %s", lines, RcFilePath);
                new BashCommand(command).logLevel(null).run();
            }
            else if(!mQuietInitialize)
            {
                CUP_LOGGER.warn("Skipping adding pyenv paths to rc file: " + RcFilePath);
            }
        }
        catch(Exception e)
        {
            CUP_LOGGER.error("Failed to update rc file({}): {}", RcFilePath, e);
            System.exit(1);
        }
    }

    @VisibleForTesting
    void installPython()
    {
        File pythonPathInPyenv = new File(String.format("%s/versions/%s/bin/python", PYENV_DIR, mPythonVersion));

        if(!pythonPathInPyenv.exists())
        {
            String command = String.format("%s install %s", PYENV_PATH, mPythonVersion);
            CUP_LOGGER.info("Installing python version {} with command: {}", mPythonVersion, command);
            new BashCommand(command).logLevel(Level.DEBUG).run();
        }
        else if(!mQuietInitialize)
        {
            CUP_LOGGER.warn("Skipping installing python {} as it exists at: {}", mPythonVersion, pythonPathInPyenv);
        }
    }

    @VisibleForTesting
    void createVirtualEnvironment()
    {
        if(!pythonPath().exists())
        {
            String command = String.format("%s virtualenv %s %s", PYENV_PATH, mPythonVersion, mVirtualEnvName);
            CUP_LOGGER.info("Creating python virtual environment with command: " + command);
            new BashCommand(command).logLevel(Level.DEBUG).run();
        }
        else if(!mQuietInitialize)
        {
            CUP_LOGGER.warn("Skipping creating python virtual environment as binary exists at: " + pythonPath());
        }
    }

    @VisibleForTesting
    void removeExistingPyenv()
    {
        CUP_LOGGER.info("Removing existing pyenv at: " + PYENV_DIR);
        try {
            FileUtils.deleteDirectory(PYENV_DIR);
        } catch(IOException e) {
            CUP_LOGGER.warn("Failed to remove pyenv: " + e);
        }
    }

    public PythonEnvironment initialize()
    {
        installPyenv();
        addPyenvPathsToRcFile();
        installPython();
        createVirtualEnvironment();
        return this;
    }

    public void pipUpgrade()
    {
        String command = "python -m pip install --upgrade pip";
        CUP_LOGGER.info("Upgrading pip with command: " + command);
        new PythonEnvCommand(this, command).showCommand().logLevel(Level.DEBUG).run();
    }

    public void pipInstall(String args)
    {
        String command = "python -m pip install " + args;
        CUP_LOGGER.info("Installing packages with command: " + command);
        new PythonEnvCommand(this, command).showCommand().logLevel(Level.DEBUG).run();
    }
}


/* setup_python_v3_9_0.sh
#!/usr/bin/env bash

source message_functions || exit 1

VER="3.9.0"

info "Installing pyenv..."
curl -L https://raw.githubusercontent.com/yyuu/pyenv-installer/master/bin/pyenv-installer | bash
ls -l ~/.pyenv/

info "Initialising pyenv..."
~/.pyenv/bin/pyenv install ${VER}
~/.pyenv/versions/${VER}/bin/pip install --upgrade pip

info "Switching to new ${VER} version"
export PATH="$HOME/.pyenv/versions/${VER}/bin:$PATH"
 */


/* do_run_cuppa.sh
#!/usr/bin/env bash

source metadata_functions || exit 1
source message_functions || exit 1
source locate_files || exit 1
source io_functions || exit 1

run_dir=$1 && shift
cuppa_jar=$1 && shift
cuppa_output_dir=$1 && shift
cuppa_mode=$1 && shift

if [[ -z "${run_dir}" || -z "${cuppa_jar}" || -z "${cuppa_output_dir}" || -z "${cuppa_mode}" ]]; then
    error "Missing params. Exiting"
fi

if [[ ${cuppa_mode} == "ALL" && ! -d "${run_dir}/isofox" ]]; then
    error "Cannot run cuppa mode ALL in case there is no isofox data! Please run dna-only"
fi

tumor_sample=$(load_tumor_sample_from_metadata ${run_dir})

tmp_sample_data_dir="${run_dir}/cuppa_sample_data_tmp"

info "Creating and populating sample data dir '${tmp_sample_data_dir}'"
create_or_cleanup_dir ${tmp_sample_data_dir}
cp -r ${run_dir}/purple/* "${tmp_sample_data_dir}/"
cp -r ${run_dir}/linx/* "${tmp_sample_data_dir}/"
cp -r ${run_dir}/virus_interpreter/* "${tmp_sample_data_dir}/"

if [[ -d "${run_dir}/isofox" ]]; then
    cp -r ${run_dir}/isofox/* "${tmp_sample_data_dir}/"
fi

create_or_cleanup_dir "${cuppa_output_dir}"

info "Running CuppaDataPrep with mode ${cuppa_mode} on ${run_dir}"
java -Xmx4G -cp ${cuppa_jar} com.hartwig.hmftools.cup.prep.CuppaDataPrep \
    -sample "${tumor_sample}" \
    -categories "${cuppa_mode}" \
    -ref_genome_version "V37" \
    -ref_alt_sj_sites "$(locate_cuppa_ref_alt_sj_sites)" \
    -sample_data_dir "${tmp_sample_data_dir}" \
    -output_dir "${cuppa_output_dir}" \
    "$@"

info "Removing temporary data dir '${tmp_sample_data_dir}'"
rm -r ${tmp_sample_data_dir}

python_version=$(python3 --version)
if [[ "${python_version}" != "Python 3.9.0" ]]; then
    info "Adding python 3.9.0 to PATH"
    export PATH="$HOME/.pyenv/versions/3.9.0/bin:$PATH"

    python_version=$(python3 --version)
    if [[ "${python_version}" != "Python 3.9.0" ]]; then
        error "Could not locate Python 3.9.0! Please run 'setup_python_v3_9_0' first. Exiting."
    fi
fi

info "Setting up pycuppa"
pycuppa_base_dir="${HOME}"
pycuppa_dir="${pycuppa_base_dir}/pycuppa"
create_or_cleanup_dir "${pycuppa_dir}"
unzip ${cuppa_jar} pycuppa/* -d "${pycuppa_base_dir}"

info "Setting up pycuppa venv"
pycuppa_venv="${pycuppa_base_dir}/pycuppa_venv"
create_or_cleanup_dir ${pycuppa_venv}
python3 -m venv ${pycuppa_venv}
source ${pycuppa_venv}/bin/activate

info "Installing pycuppa tool into venv"
pip3 install --upgrade pip
pip3 install "${pycuppa_dir}"

info "Running pycuppa DNA on ${run_dir}"
python3 -m cuppa.predict \
    --sample_id ${tumor_sample} \
    --features_path "$(locate_cuppa_features_tsv ${cuppa_output_dir})" \
    --classifier_path "$(locate_cuppa_classifier_pickle)" \
    --cv_predictions_path "$(locate_cuppa_cv_predictions)" \
    --output_dir "${cuppa_output_dir}"

info "Regenerating CUPPA R visualization for ${sample} in ${cuppa_output_dir} using pilot visualizer"
Rscript --vanilla $(locate_pilot_cuppa_visualizer) \
    "$(locate_cuppa_vis_data_tsv ${cuppa_output_dir})" \
    "${cuppa_output_dir}/${tumor_sample}.cuppa.pilot.vis.png"
 */
