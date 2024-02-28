package com.hartwig.hmftools.esvee.alignment;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import com.hartwig.hmftools.esvee.SvConfig;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;

public class BwaAligner implements Aligner
{
    private final BwaMemIndex mIndex;
    private final BwaMemAligner mAligner;

    public BwaAligner(final SvConfig config)
    {
        if(!config.RefGenomeImageFile.isEmpty() && Files.exists(Paths.get(config.RefGenomeImageFile)))
        {
            loadAlignerLibrary();

            mIndex = new BwaMemIndex(config.RefGenomeImageFile);
            mAligner = new BwaMemAligner(mIndex);
        }
        else
        {
            mIndex = null;
            mAligner = null;
        }
    }

    private void loadAlignerLibrary()
    {
        final var props = System.getProperties();
        final String candidateBWAPath = "libbwa." + props.getProperty("os.arch") + osExtension();

        if(System.getProperty("LIBBWA_PATH") == null && new File(candidateBWAPath).exists())
            System.setProperty("LIBBWA_PATH", new File(candidateBWAPath).getAbsolutePath());
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
