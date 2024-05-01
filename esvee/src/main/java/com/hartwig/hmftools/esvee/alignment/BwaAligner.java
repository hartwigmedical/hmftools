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

    public static void loadAlignerLibrary(@Nullable final String bwaLibraryPath)
    {
        final var props = System.getProperties();

        if(bwaLibraryPath != null)
        {
            System.setProperty("LIBBWA_PATH", new File(bwaLibraryPath).getAbsolutePath());
            return;
        }

        String osName = System.getProperty("os.name");
        String osLibExtension = osExtension(osName);
        String osArchitecture = System.getProperty("os.arch");

        String candidateBWAPath = null;

        if(osName.contains("Mac") && osArchitecture.equals("aarch64"))
        {
            // candidateBWAPath = Resources.getResource("libbwa.Darwin.dylib").getPath();
        }
        else
        {
            candidateBWAPath = "libbwa." + props.getProperty("os.arch") + osLibExtension;
        }

        if(System.getProperty("LIBBWA_PATH") == null && candidateBWAPath != null && new File(candidateBWAPath).exists())
        {
            System.setProperty("LIBBWA_PATH", new File(candidateBWAPath).getAbsolutePath());
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
