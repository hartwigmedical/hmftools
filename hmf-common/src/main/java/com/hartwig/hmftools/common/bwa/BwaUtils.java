package com.hartwig.hmftools.common.bwa;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;

import org.jetbrains.annotations.Nullable;

public final class BwaUtils
{
    public static final String LIBBWA_PATH = "LIBBWA_PATH"; // as expected by the BWA library
    public static final String LIBBWA_PREFIX = "libbwa.";

    public static final String MAC_OS = "Mac";
    public static final String MAC_ARCH = "aarch64";
    public static final String MAC_BWA_LIB = "libbwa.Darwin.dylib";

    public static final String BWA_LIB_PATH = "bwa_lib";
    public static final String BWA_LIB_PATH_DESC = "Path to BWA library";

    public static void loadAlignerLibrary(@Nullable final String bwaLibraryPath)
    {
        if(System.getProperty(LIBBWA_PATH) != null)
            return;

        if(bwaLibraryPath != null)
        {
            System.setProperty(LIBBWA_PATH, new File(bwaLibraryPath).getAbsolutePath());
            return;
        }

        String osName = System.getProperty("os.name");
        String osLibExtension = osExtension(osName);
        String osArchitecture = System.getProperty("os.arch");

        String candidateBWAPath = null;

        if(osName.contains(MAC_OS) && osArchitecture.equals(MAC_ARCH))
        {
            candidateBWAPath = MAC_BWA_LIB;
        }
        else
        {
            candidateBWAPath = LIBBWA_PREFIX + osArchitecture + osLibExtension;
        }

        if(Files.exists(Paths.get(candidateBWAPath)))
        {
            System.setProperty(LIBBWA_PATH, new File(candidateBWAPath).getAbsolutePath());
        }
    }

    private static String osExtension(final String osName)
    {
        if(osName.contains("Mac"))
            return ".dylib";
        else if(osName.contains("Win"))
            return ".dll";
        else
            return ".so";
    }
}
