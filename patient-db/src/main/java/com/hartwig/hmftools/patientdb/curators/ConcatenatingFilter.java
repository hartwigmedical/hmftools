package com.hartwig.hmftools.patientdb.curators;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;

import org.apache.lucene.analysis.TokenFilter;
import org.apache.lucene.analysis.TokenStream;
import org.apache.lucene.analysis.tokenattributes.CharTermAttribute;
import org.apache.lucene.analysis.tokenattributes.OffsetAttribute;
import org.apache.lucene.analysis.tokenattributes.PositionIncrementAttribute;
import org.apache.lucene.analysis.tokenattributes.PositionLengthAttribute;
import org.apache.lucene.analysis.tokenattributes.TypeAttribute;
import org.apache.lucene.util.AttributeSource;

public class ConcatenatingFilter extends TokenFilter {
    private final CharTermAttribute termAttribute = addAttribute(CharTermAttribute.class);
    private final OffsetAttribute offsetAtt = addAttribute(OffsetAttribute.class);
    private final PositionIncrementAttribute posIncrAtt = addAttribute(PositionIncrementAttribute.class);
    private final PositionLengthAttribute posLenAtt = addAttribute(PositionLengthAttribute.class);
    private final TypeAttribute typeAtt = addAttribute(TypeAttribute.class);

    private List<String> terms = null;
    private AttributeSource.State finalState;

    private final char separator;
    private boolean inputEnded = false;

    ConcatenatingFilter(TokenStream input, char separator) {
        super(input);
        this.separator = separator;
    }

    @Override
    public final boolean incrementToken() throws IOException {
        if (terms != null) {
            return false;
        }
        boolean result = buildSingleOutputToken();
        finalState = captureState();
        return result;
    }

    private boolean buildSingleOutputToken() throws IOException {
        inputEnded = false;
        terms = Lists.newArrayList();

        while (input.incrementToken()) {
            terms.add(termAttribute.toString());
        }

        input.end();
        inputEnded = true;

        //MIVO: set attributes for the single output token
        offsetAtt.setOffset(0, offsetAtt.endOffset());
        posLenAtt.setPositionLength(1);
        posIncrAtt.setPositionIncrement(1);
        typeAtt.setType("concatenation");

        if (terms.size() < 1) {
            termAttribute.setEmpty();
            return false;
        }

        StringBuilder sb = new StringBuilder();
        for (String term : terms) {
            if (sb.length() >= 1) {
                sb.append(separator);
            }
            sb.append(term);
        }

        termAttribute.setEmpty().append(sb);
        terms = Lists.newArrayList();
        return true;
    }

    @Override
    public final void end() throws IOException {
        if (!inputEnded) {
            // Rare case - If an IOException occurs while performing buildSingleOutputToken
            // we may not have called input.end() already
            input.end();
            inputEnded = true;
        }

        if (finalState != null) {
            restoreState(finalState);
        }
    }

    @Override
    public void reset() throws IOException {
        super.reset();
        inputEnded = false;
        terms = null;
    }
}