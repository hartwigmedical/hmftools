package com.hartwig.hmftools.esvee.common.saga;

import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;

import java.util.List;

public record SagaAssembly(
        String fastaLabel,
        SagaVariant variant,
        // Junctions occur just before each of these indices in the assembly sequence.
        // Usually length 2. Can be length 1 for DELs where the junction is a single position.
        // E.g.
        // sequence = RRRJJJRRR
        // junctionOffsets[0] = 3
        // junctionOffsets[1] = 6
        List<Integer> junctionOffsets,
        String sequence
)
{
    public SagaAssembly
    {
        if(junctionOffsets.size() != 1 && junctionOffsets.size() != 2)
        {
            throw new IllegalArgumentException("Expected 1 or 2 junction offsets");
        }

        if(variant.isInsert() && junctionOffsets.size() != 2)
        {
            throw new IllegalArgumentException("Insert should have 2 junction offsets");
        }

        if(!junctionOffsets.stream().allMatch(offset -> offset >= 1 && offset < sequence.length()))
        {
            throw new IllegalArgumentException("Junction offsets out of bounds");
        }

        for(int i = 0; i < junctionOffsets.size() - 1; i++)
        {
            if(junctionOffsets.get(i) >= junctionOffsets.get(i + 1))
            {
                throw new IllegalArgumentException("Junction offsets not in ascending order");
            }
        }
    }

    public int length()
    {
        return sequence.length();
    }

    public List<SagaJunctionInfo> junctions()
    {
        return junctionOffsets.stream().map(SagaJunctionInfo::new).toList();
    }

    public static SagaAssembly fromFastaRecord(final String label, final String sequence)
    {
        // E.g.:
        // SvimAsm00000237|chr1:181626:1|chr1:181627:-1|150|285
        // SvimAsm00000238|chr1:368909:1|chr1:369380:-1|150|

        String[] parts = label.split("\\|");
        if(!(parts.length == 4 || parts.length == 5))
        {
            SV_LOGGER.error("Expected 4 or 5 parts but got {}", parts.length);
            throw new IllegalArgumentException("Invalid fasta label: " + label);
        }
        String id = parts[0];
        SagaBreakend breakend1 = SagaBreakend.fromString(parts[1]);
        SagaBreakend breakend2 = SagaBreakend.fromString(parts[2]);
        int junctionOffset1 = Integer.parseInt(parts[3]);
        Integer junctionOffset2 = parts.length >= 5 ? Integer.parseInt(parts[4]) : null;
        if(junctionOffset2 != null && junctionOffset1 >= junctionOffset2)
        {
            throw new IllegalArgumentException("Invalid junction offsets");
        }
        String insertSequence = sequence.substring(junctionOffset1, junctionOffset2 == null ? junctionOffset1 : junctionOffset2);
        SagaVariant variant = new SagaVariant(id, breakend1, breakend2, insertSequence);
        List<Integer> junctionOffsets = junctionOffset2 == null ? List.of(junctionOffset1) : List.of(junctionOffset1, junctionOffset2);
        return new SagaAssembly(label, variant, junctionOffsets, sequence);
    }
}
