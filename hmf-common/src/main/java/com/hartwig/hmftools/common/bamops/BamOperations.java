package com.hartwig.hmftools.common.bamops;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public final class BamOperations
{
    private static final String INDEX_COMMAND = "index";
    private static final String MERGE_COMMAND = "merge";
    private static final String CONCATENATE_COMMAND = "cat";
    private static final String SORT_COMMAND = "sort";

    protected static final Logger BOP_LOGGER = LogManager.getLogger(BamOperations.class);

    public static boolean mergeBams(
            final BamToolName toolName, final String toolPath, final String outputBam, final List<String> inputBams, final int threads)
    {
        BOP_LOGGER.trace("merging {} bams", inputBams.size());

        List<String> commandArgs = Lists.newArrayList();

        commandArgs.add(toolPath);
        commandArgs.add(MERGE_COMMAND);

        addThreadsArg(toolName, commandArgs, threads);
        commandArgs.add(outputBam);

        for(String threadBam : inputBams)
        {
            commandArgs.add(threadBam);
        }

        if(!executeCommand(commandArgs, outputBam))
            return false;

        BOP_LOGGER.trace("merge complete");
        return true;
    }

    public static boolean concatenateBams(
            final BamToolName toolName, final String toolPath, final String outputBam, final List<String> inputBams, final int threads)
    {
        if(toolName != BamToolName.SAMTOOLS)
        {
            BOP_LOGGER.error("{} concatenation no supported", toolName);
            return false;
        }

        BOP_LOGGER.trace("concatenating {} bams", inputBams.size());

        List<String> commandArgs = Lists.newArrayList();

        commandArgs.add(toolPath);
        commandArgs.add(CONCATENATE_COMMAND);

        addThreadsArg(toolName, commandArgs, threads);

        commandArgs.add("-o");
        commandArgs.add(outputBam);

        commandArgs.add("--no-PG");

        for(String threadBam : inputBams)
        {
            commandArgs.add(threadBam);
        }

        if(!executeCommand(commandArgs, outputBam))
            return false;

        return true;
    }

    public static boolean indexBam(final BamToolName toolName, final String toolPath, final String bamFilename, final int threads)
    {
        BOP_LOGGER.trace("indexing bam({})", bamFilename);

        List<String> commandArgs = Lists.newArrayList();

        commandArgs.add(toolPath);
        commandArgs.add(INDEX_COMMAND);
        addThreadsArg(toolName, commandArgs, threads);
        commandArgs.add(bamFilename);

        if(!executeCommand(commandArgs, bamFilename))
            return false;

        BOP_LOGGER.trace("index complete");
        return true;
    }

    public static boolean sortBam(
            final BamToolName toolName, final String toolPath, final String inputBam, final String outputBam, final int threads)
    {
        List<String> commandArgs = Lists.newArrayList();

        commandArgs.add(toolPath);
        commandArgs.add(SORT_COMMAND);
        addThreadsArg(toolName, commandArgs, threads);

        // default memory per thread according to samtools doco is 768MB, could configure as a function of max heap used by MarkDups
        // commandArgs.add("-m");
        // commandArgs.add("1G");

        // commandArgs.add("-T");
        // commandArgs.add("tmp"); // these are default

        if(toolName == BamToolName.SAMTOOLS)
        {
            commandArgs.add("-O");
            commandArgs.add("bam");
        }

        commandArgs.add(inputBam);

        commandArgs.add("-o");
        commandArgs.add(outputBam);

        if(!executeCommand(commandArgs, outputBam))
            return false;

        BOP_LOGGER.trace("sort complete");
        return true;
    }


    private static void addThreadsArg(final BamToolName toolName, final List<String> commandArgs, final int threads)
    {
        if(threads > 1)
        {
            commandArgs.add(toolName.toolThreadArgument());
            commandArgs.add(String.valueOf(threads));
        }
    }

    private static boolean executeCommand(final List<String> commandArgs, final String outputPrefix)
    {
        String redirectOutputFile = outputPrefix + ".out";
        String redirectErrFile = outputPrefix + ".err";

        String[] command = new String[commandArgs.size()];
        for(int i = 0; i < commandArgs.size(); ++i)
        {
            command[i] = commandArgs.get(i);
        }

        try
        {
            Process process = new ProcessBuilder(command)
                    .redirectOutput(new File(redirectOutputFile))
                    .redirectError(new File(redirectErrFile))
                    .start();

            int result = process.waitFor();

            if(result != 0)
            {
                BOP_LOGGER.error("error running command({}) for file({})", commandToStr(command), outputPrefix);
                logErrorFile(redirectErrFile);
                // System.err.print(new String(process.getErrorStream().readAllBytes()));
                return false;
            }

            // clean-up process log files
            Files.deleteIfExists(Paths.get(redirectOutputFile));
            Files.deleteIfExists(Paths.get(redirectErrFile));
        }
        catch(Exception e)
        {
            BOP_LOGGER.error("error running command({}) for file({}): {}", commandToStr(command), outputPrefix, e.toString());
            return false;
        }

        return true;
    }

    private static String commandToStr(final String[] command)
    {
        return Arrays.stream(command).collect(Collectors.joining(" "));
    }

    private static void logErrorFile(final String errorFile)
    {
        try
        {
            BOP_LOGGER.error("error file({}) contents", errorFile);

            for(String line : Files.readAllLines(Paths.get(errorFile)))
            {
                BOP_LOGGER.error("{}", line);
            }
        }
        catch(Exception e)
        {
            BOP_LOGGER.error("cannot read  error file: {}", e.toString());
        }
    }
}
