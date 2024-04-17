package com.hartwig.hmftools.esvee.alignment;

import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.List;

import com.hartwig.hmftools.esvee.AssemblyConfig;

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

            // TEMP: until can resolve local ARM library issues
            try
            {
                index = new BwaMemIndex(refGenomeImageFile);
            }
            catch(Exception e)
            {
                SV_LOGGER.error("failed to initialise BWA aligner: {}", e.toString());
            }

            mAligner = index != null ? new BwaMemAligner(index) : null;
        }
        else
        {
            mAligner = null;
        }
    }

    public static void loadAlignerLibrary(@Nullable final String bwaLibPath)
    {
        final var props = System.getProperties();
        String candidateBWAPath = bwaLibPath != null ? bwaLibPath : "libbwa." + props.getProperty("os.arch") + osExtension();

        if(System.getProperty("LIBBWA_PATH") == null && new File(candidateBWAPath).exists())
        {
            System.setProperty("LIBBWA_PATH", new File(candidateBWAPath).getAbsolutePath());
        }
    }

    private static String osExtension()
    {
        final String osName = System.getProperty("os.name");
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

        /*
        @Nullable
        final BwaMemAlignment best = alignmentSet.stream()
                .filter(alignment -> alignment.getRefStart() != -1)
                .max(Comparator.comparingInt(BwaMemAlignment::getMapQual)
                        .thenComparingInt(BwaMemAlignment::getAlignerScore)
                        .thenComparing(Comparator.comparingInt(BwaMemAlignment::getNMismatches).reversed()))
                .orElse(null);
        if(best == null)
            return List.of(com.hartwig.hmftools.esvee.old.Alignment.unmapped(sequence.getLength()));

        if(SAMFlag.getFlags(best.getSamFlag()).contains(SAMFlag.READ_UNMAPPED))
            return List.of(com.hartwig.hmftools.esvee.old.Alignment.unmapped(sequence.getLength()));
        */
    }

}
