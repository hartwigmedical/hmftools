package com.hartwig.hmftools.esvee.alignment;

import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_INDEL_LENGTH;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.List;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;
import org.jetbrains.annotations.Nullable;

public class BwaAligner implements Aligner
{
    private final BwaMemAligner mAligner;

    public BwaAligner(final String refGenomeImageFile)
    {
        if(!refGenomeImageFile.isEmpty() && Files.exists(Paths.get(refGenomeImageFile)))
        {
            BwaMemIndex index = null;

            try
            {
                index = new BwaMemIndex(refGenomeImageFile);
            }
            catch(Exception e)
            {
                SV_LOGGER.error("failed to initialise BWA aligner: {}", e.toString());
            }

            if(index != null)
            {
                mAligner = new BwaMemAligner(index);
                mAligner.setBandwidthOption(MIN_INDEL_LENGTH - 1);
            }
            else
            {
                mAligner = null;
            }
        }
        else
        {
            mAligner = null;
        }
    }

    private static final String LIBBWA_PATH = "LIBBWA_PATH"; // as expected by the BWA library
    private static final String LIBBWA_PREFIX = "libbwa.";

    private static final String MAC_OS = "Mac";
    private static final String MAC_ARCH = "aarch64";
    private static final String MAC_BWA_LIB = "libbwa.Darwin.dylib";

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

    @Override
    public List<BwaMemAlignment> alignSequence(final byte[] bases)
    {
        if(mAligner == null)
            return Collections.emptyList();

        List<BwaMemAlignment> alignmentSet = mAligner.alignSeqs(List.of(bases)).get(0);

        return alignmentSet;
    }
}
