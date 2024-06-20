package com.hartwig.hmftools.common.utils.r;

import static java.lang.String.format;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.Arrays;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.IOUtil;

/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

public final class RExecutor
{
    private static final String R_EXE = "Rscript";
    private static final Logger LOGGER = LogManager.getLogger(RExecutor.class);

    public static int executeFromClasspath(final String rScriptName, final String... arguments)
            throws IOException, InterruptedException
    {
        final File scriptFile = writeScriptFile(rScriptName);

        final int returnCode = executeFromFile(rScriptName, scriptFile, arguments);
        htsjdk.samtools.util.IOUtil.deleteFiles(scriptFile);
        return returnCode;
    }

    private static int executeFromFile(final String rScriptName, final File scriptFile, final String... arguments)
            throws IOException, InterruptedException
    {
        final String[] command = new String[arguments.length + 2];
        command[0] = R_EXE;
        command[1] = scriptFile.getAbsolutePath();
        System.arraycopy(arguments, 0, command, 2, arguments.length);

        final File outputFile = File.createTempFile(rScriptName, ".out");

        LOGGER.debug(format("executing R script via command: %s", CollectionUtil.join(Arrays.asList(command), " ")));
        Process process = new ProcessBuilder(command).redirectOutput(outputFile).start();

        int result = process.waitFor();

        if(result != 0)
        {
            System.err.print(new String(process.getErrorStream().readAllBytes()));
            LOGGER.fatal("error executing R script");
        }

        return result;
    }

    private static File writeScriptFile(final String rScriptName) throws IOException
    {
        InputStream scriptStream = null;
        OutputStream scriptFileStream = null;
        try
        {
            scriptStream = RExecutor.class.getClassLoader().getResourceAsStream(rScriptName);
            if(scriptStream == null)
            {
                throw new IllegalArgumentException("Script [" + rScriptName + "] not found in classpath");
            }
            final File scriptFile = File.createTempFile("script", ".R");
            scriptFileStream = IOUtil.openFileForWriting(scriptFile);
            IOUtil.copyStream(scriptStream, scriptFileStream);
            return scriptFile;
        }
        finally
        {
            if(scriptStream != null)
            {
                try
                {
                    scriptStream.close();
                }
                catch(IOException ignored)
                {
                }
            }
            if(scriptFileStream != null)
            {
                try
                {
                    scriptFileStream.close();
                }
                catch(IOException ignored)
                {
                }
            }
        }
    }
}
