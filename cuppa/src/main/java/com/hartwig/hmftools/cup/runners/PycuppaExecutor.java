package com.hartwig.hmftools.cup.runners;

import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.io.Resources;
import org.apache.commons.io.FileUtils;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;

import jdk.jfr.Experimental;

@Experimental
public class PycuppaExecutor
{
    private static final String[] SHELL_PATHS_TO_CHECK = {"/bin/bash", "/bin/sh", "/bin/zsh"};

    private static final File TMP_DIR = new File(System.getProperty("java.io.tmpdir"));
    private static final File PYCUPPA_TMP_DIR = new File(TMP_DIR + "/pycuppa/");
    private static final File PYCUPPA_RESOURCE_DIR = new File(Resources.getResource("pycuppa/").getPath());

    public final File mVirtualEnvDir;
    private final File mVirtualEnvActivator;

    public PycuppaExecutor(String virtualEnvDir)
    {
        mVirtualEnvDir = new File(virtualEnvDir);
        mVirtualEnvActivator = new File(mVirtualEnvDir + "/bin/activate");
    }

    private static String getShellPath() throws FileNotFoundException
    {
        String shellPath = null;
        for(String path : SHELL_PATHS_TO_CHECK)
        {
            if(new File(path).exists())
            {
                shellPath = path;
                break;
            }
        }

        if(shellPath == null)
            throw new FileNotFoundException("Could not find shell path. Tried the paths: " + Arrays.toString(SHELL_PATHS_TO_CHECK));

        return shellPath;
    }

    public static List<String> runBashCommand(String command, Level logLevel)
    {
        List<String> stdout = new ArrayList<>();
        int exitCode = 0;
        try
        {
            ProcessBuilder processBuilder = new ProcessBuilder(getShellPath(), "-c", command);
            processBuilder.redirectErrorStream(true); // Merge stdout and stderr

            Process process = processBuilder.start();

            BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
            String line;
            while( (line = reader.readLine()) != null)
            {
                stdout.add(line);

                if(logLevel != Level.OFF) // Added if clause because logging still prints at Level.OFF
                    CUP_LOGGER.log(logLevel, line);
            }

            exitCode = process.waitFor();
            if(exitCode > 0)
                throw new RuntimeException();
        }
        catch(Exception e)
        {
            CUP_LOGGER.error("Bash command({}) failed with exit code {}: {}", command, String.valueOf(exitCode), e);
            System.exit(exitCode);
        }

        return stdout;
    }

    public List<String> runBashCommandInVirtualEnv(String command, Level logLevel)
    {
        return runBashCommand(String.format("source %s && %s && deactivate", mVirtualEnvActivator, command), logLevel);
    }

    private void createVirtualEnvIfMissing()
    {
        if(mVirtualEnvDir.exists())
        {
            if(!mVirtualEnvActivator.exists())
            {
                CUP_LOGGER.error("dir({}) exists but is not a python virtual environment", mVirtualEnvDir);
                System.exit(1);
            }

            // CUP_LOGGER.debug("Using existing python virtual environment: " + mVirtualEnvDir);
            return;
        }

        runBashCommand("python3 -m venv " + mVirtualEnvDir, Level.DEBUG);
        CUP_LOGGER.info("Created python virtual environment at: " + mVirtualEnvDir);
    }

    private static void extractPycuppaToTmpDir() throws IOException, URISyntaxException
    {
        if(PYCUPPA_TMP_DIR.exists())
        {
            CUP_LOGGER.debug("Removing existing pycuppa tmp dir({})", PYCUPPA_TMP_DIR);
            FileUtils.deleteDirectory(PYCUPPA_TMP_DIR);
        }

        if(PYCUPPA_RESOURCE_DIR.exists())
        {
            // This clause allows PycuppaExecutor to work when run from intellij
            CUP_LOGGER.debug(" Copying pycuppa from resource path({}) to tmp dir({})", PYCUPPA_RESOURCE_DIR, PYCUPPA_TMP_DIR);
            FileUtils.copyDirectory(PYCUPPA_RESOURCE_DIR, PYCUPPA_TMP_DIR);
        }
        else
        {
            String cuppaJarPath = new File( PycuppaExecutor.class.getProtectionDomain().getCodeSource().getLocation().toURI() ).getPath();

            CUP_LOGGER.debug(" Extracting pycuppa from jar({}) to tmp dir({})", cuppaJarPath, PYCUPPA_TMP_DIR);

            // Unzip with `jar` rather than `unzip` to avoid unzipping the whole jar
            String extractCommand = String.format("cd %s && jar -xf %s pycuppa", TMP_DIR, cuppaJarPath);
            runBashCommand(extractCommand, Level.DEBUG);
        }

        if(!PYCUPPA_TMP_DIR.exists())
            throw new RuntimeException("Failed to extract pycuppa to tmp dir(" + PYCUPPA_TMP_DIR + ")");
    }

    private boolean pycuppaInstalled()
    {
        List<String> pipOutput = runBashCommandInVirtualEnv("pip --disable-pip-version-check list", Level.OFF);

        boolean pycuppaExists = false;
        for(String line : pipOutput)
        {
            if(line.startsWith("pycuppa"))
            {
                pycuppaExists = true;
                break;
            }
        }

        return pycuppaExists;
    }

    private void installPycuppaIfMissing()
    {
        if(pycuppaInstalled())
            return;

        try
        {
            CUP_LOGGER.info("Installing pycuppa to virtual environment: " + mVirtualEnvDir);

            extractPycuppaToTmpDir();
            runBashCommandInVirtualEnv("pip --disable-pip-version-check install " + PYCUPPA_TMP_DIR, Level.DEBUG);
        }
        catch(Exception e)
        {
            CUP_LOGGER.error("Failed to install pycuppa: ", e);
        }
        finally
        {
            try
            {
                FileUtils.deleteDirectory(PYCUPPA_TMP_DIR);
            }
            catch(IOException e)
            {
                CUP_LOGGER.error("Failed to delete tmp pycuppa dir({})", PYCUPPA_TMP_DIR);
                System.exit(1);
            }
        }
    }

    public void initialize()
    {
        createVirtualEnvIfMissing();
        installPycuppaIfMissing();
    }

    @VisibleForTesting
    public static void main(final String[] args) throws Exception
    {
        Configurator.setLevel(CUP_LOGGER.getName(), Level.DEBUG);
        PycuppaExecutor executor = new PycuppaExecutor(args[0]);
        executor.initialize();
    }
}
