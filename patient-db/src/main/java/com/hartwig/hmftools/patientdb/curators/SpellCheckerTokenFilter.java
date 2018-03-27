package com.hartwig.hmftools.patientdb.curators;

import java.io.IOException;

import org.apache.lucene.analysis.TokenFilter;
import org.apache.lucene.analysis.TokenStream;
import org.apache.lucene.analysis.tokenattributes.CharTermAttribute;
import org.apache.lucene.search.spell.SpellChecker;
import org.jetbrains.annotations.NotNull;

class SpellCheckerTokenFilter extends TokenFilter {
    private final SpellChecker spellChecker;
    private final CharTermAttribute charTermAttribute = addAttribute(CharTermAttribute.class);

    SpellCheckerTokenFilter(@NotNull final TokenStream tokenStream, @NotNull final SpellChecker spellChecker) {
        super(tokenStream);
        this.spellChecker = spellChecker;
    }

    @Override
    public final boolean incrementToken() throws IOException {
        if (!this.input.incrementToken()) {
            return false;
        }
        final String currentTokenInStream = this.input.getAttribute(CharTermAttribute.class).toString();
        final String nextToken = spellcheck(currentTokenInStream);
        this.charTermAttribute.setEmpty().append(nextToken);
        return true;
    }

    @NotNull
    private String spellcheck(@NotNull final String searchString) throws IOException {
        if (spellChecker.exist(searchString)) {
            return searchString;
        } else {
            String[] suggestions = spellChecker.suggestSimilar(searchString, 20);
            if (suggestions.length > 0) {
                return suggestions[0];
            } else {
                return searchString;
            }
        }
    }
}
