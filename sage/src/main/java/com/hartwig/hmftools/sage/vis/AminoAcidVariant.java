package com.hartwig.hmftools.sage.vis;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.codon.AminoAcids.TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

public interface AminoAcidVariant
{
    char INDEL = '.';

    // TODO: do we need to parse this?
    record AminoAcidSubstitution(int aminoAcidPos, char ref, char alt) implements AminoAcidVariant
    {
        private static final Pattern PATTERN = Pattern.compile("^p\\.([^0-9]+)([0-9]+)([^0-9]+)$");

        public static AminoAcidSubstitution parse(final String s)
        {
            Matcher matcher = PATTERN.matcher(s);
            if(!matcher.find())
                throw new RuntimeException(format("%s doesn't match the expected format for a HGVS", s));

            String triRef = matcher.group(1);
            int pos = Integer.parseInt(matcher.group(2));
            String triAlt = matcher.group(3);

            String ref = TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.get(triRef);
            if(triAlt.equals("del"))
                return new AminoAcidSubstitution(pos, ref.charAt(0), INDEL);

            String alt = TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.get(triAlt);
            return new AminoAcidSubstitution(pos, ref.charAt(0), alt.charAt(0));
        }
    }
}
