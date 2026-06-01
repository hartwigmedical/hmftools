package com.hartwig.hmftools.esvee.common;

import static java.lang.String.format;

import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.BasePosition;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class SagaResource
{
    private final String mFastaPath;
    private final List<AssemblyMetadata> mAssemblies;
    private final SAMSequenceDictionary mSequenceDictionary;

    public SagaResource(final String fastaPath)
    {
        mFastaPath = fastaPath;
        mAssemblies = new ArrayList<>(100_000);
        mSequenceDictionary = loadAssemblies(fastaPath, mAssemblies);
    }

    public String bwaIndexImagePath()
    {
        return mFastaPath + ".img";
    }

    public SAMSequenceDictionary samDict()
    {
        return mSequenceDictionary;
    }

    public List<AssemblyMetadata> assemblies()
    {
        return mAssemblies;
    }

    public record Breakend(
            BasePosition position,
            Orientation orientation
    )
    {
        public static Breakend fromString(final String string)
        {
            String[] parts = string.split(":");
            if(parts.length != 3)
            {
                throw new IllegalArgumentException("Invalid breakend string: " + string);
            }
            String chromosome = parts[0];
            int position = Integer.parseInt(parts[1]);
            Orientation orientation = Orientation.fromByteStr(parts[2]);
            return new Breakend(new BasePosition(chromosome, position), orientation);
        }

        @Override
        @NotNull
        public String toString()
        {
            return format("%s:%s", position, orientation);
        }
    }

    public record Variant(
            String id,
            Breakend breakend1,
            Breakend breakend2
    )
    {
        public Stream<Breakend> breakends()
        {
            return Stream.of(breakend1, breakend2);
        }

        @Override
        @NotNull
        public String toString()
        {
            return format("%s %s-%s", id, breakend1, breakend2);
        }

        public String noSpaceString()
        {
            String result = format("%s:%s-%s", id, breakend1, breakend2);
            result = result.replace(" ", "_");
            return result;
        }
    }

    public record AssemblyMetadata(
            String fastaLabel,
            Variant variant,
            // Junctions occur just before each of these indices in the assembly sequence.
            // Usually length 2. Can be length 1 for DELs where the junction is a single position.
            List<Integer> junctionOffsets,
            int assemblyLength
    )
    {
        public AssemblyMetadata
        {
            // E.g.
            // sequence = RRRJJJRRR
            // junctionOffset[0] = 3
            // junctionOffset[1] = 6
            if(!junctionOffsets.stream().allMatch(offset -> offset >= 1 && offset < assemblyLength))
            {
                throw new IllegalArgumentException("Junction offsets out of bounds");
            }
        }

        public String variantId()
        {
            return variant.id();
        }

        public static AssemblyMetadata fromFastaLabel(final String fastaLabel, int assemblyLength)
        {
            // E.g.:
            // SvimAsm00000237|chr1:181626:1|chr1:181627:-1|150|285
            // SvimAsm00000238|chr1:368909:1|chr1:369380:-1|150|

            String[] parts = fastaLabel.split("\\|");
            if(!(parts.length == 4 || parts.length == 5))
            {
                SV_LOGGER.error("Expected 4 or 5 parts but got {}", parts.length);
                throw new IllegalArgumentException("Invalid fasta label: " + fastaLabel);
            }
            String id = parts[0];
            Breakend breakend1 = Breakend.fromString(parts[1]);
            Breakend breakend2 = Breakend.fromString(parts[2]);
            int junctionOffset1 = Integer.parseInt(parts[3]);
            Integer junctionOffset2 = parts.length >= 5 ? Integer.parseInt(parts[4]) : null;
            if(junctionOffset2 != null && junctionOffset1 >= junctionOffset2)
            {
                throw new IllegalArgumentException("Invalid junction offsets");
            }
            Variant variant = new Variant(id, breakend1, breakend2);
            List<Integer> junctionOffsets = junctionOffset2 == null ? List.of(junctionOffset1) : List.of(junctionOffset1, junctionOffset2);
            return new AssemblyMetadata(fastaLabel, variant, junctionOffsets, assemblyLength);
        }
    }

    private static SAMSequenceDictionary loadAssemblies(final String fastaPath, List<AssemblyMetadata> assemblies)
    {
        SAMSequenceDictionary samDict;
        try(IndexedFastaSequenceFile fasta = loadFasta(fastaPath))
        {
            samDict = fasta.getSequenceDictionary();
            samDict.getSequences().stream()
                    .map(seq -> AssemblyMetadata.fromFastaLabel(seq.getSequenceName(), seq.getSequenceLength()))
                    .forEach(assemblies::add);
        }
        catch(IOException e)
        {
            throw new RuntimeException("Failed to load SAGA resource", e);
        }

        // Variants are keyed by ID later, so ensure the IDs are unique.
        long uniqueVariantIdCount = assemblies.stream().map(assembly -> assembly.variant().id()).distinct().count();
        if(uniqueVariantIdCount != assemblies.size())
        {
            throw new RuntimeException("Duplicate variant IDs in SAGA resource");
        }

        SV_LOGGER.debug("loaded {} SAGA variant assemblies from file({})", assemblies.size(), fastaPath);

        return samDict;
    }

    private static IndexedFastaSequenceFile loadFasta(final String fastaPath) throws IOException
    {
        return new IndexedFastaSequenceFile(new File(fastaPath));
    }
}
