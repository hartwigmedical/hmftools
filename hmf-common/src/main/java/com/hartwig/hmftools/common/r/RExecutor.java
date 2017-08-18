package com.hartwig.hmftools.common.r;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.Arrays;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.IOUtil;

public class RExecutor {

    private static final String R_EXE = "Rscript";
    private static final Logger LOGGER = LogManager.getLogger(RExecutor.class);

    public static int executeFromClasspath(final String rScriptName, final String... arguments) throws IOException, InterruptedException {
        final File scriptFile = writeScriptFile(rScriptName);
        final int returnCode = executeFromFile(scriptFile, arguments);
        htsjdk.samtools.util.IOUtil.deleteFiles(scriptFile);
        return returnCode;
    }

    public static int executeFromFile(final File scriptFile, final String... arguments) throws IOException, InterruptedException {
        final String[] command = new String[arguments.length + 2];
        command[0] = R_EXE;
        command[1] = scriptFile.getAbsolutePath();
        System.arraycopy(arguments, 0, command, 2, arguments.length);
        LOGGER.info(String.format("Executing R script via command: %s", CollectionUtil.join(Arrays.asList(command), " ")));

        return Runtime.getRuntime().exec(command).waitFor();
    }

    private static File writeScriptFile(final String rScriptName) throws IOException {
        InputStream scriptStream = null;
        OutputStream scriptFileStream = null;
        try {
            scriptStream = RExecutor.class.getClassLoader().getResourceAsStream(rScriptName);
            if (scriptStream == null) {
                throw new IllegalArgumentException("Script [" + rScriptName + "] not found in classpath");
            }
            final File scriptFile = File.createTempFile("script", ".R");
            scriptFileStream = IOUtil.openFileForWriting(scriptFile);
            IOUtil.copyStream(scriptStream, scriptFileStream);
            return scriptFile;
        } finally {
            if (scriptStream != null) {
                try {
                    scriptStream.close();
                } catch (IOException ignored) {
                }
            }
            if (scriptFileStream != null) {
                try {
                    scriptFileStream.close();
                } catch (IOException ignored) {
                }
            }
        }
    }
}
