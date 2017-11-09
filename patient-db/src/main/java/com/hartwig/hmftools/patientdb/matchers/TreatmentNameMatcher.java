package com.hartwig.hmftools.patientdb.matchers;

import java.io.File;
import java.io.IOException;
import java.io.StringReader;
import java.nio.charset.Charset;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.apache.lucene.analysis.Analyzer;
import org.apache.lucene.analysis.LowerCaseFilter;
import org.apache.lucene.analysis.TokenFilter;
import org.apache.lucene.analysis.TokenStream;
import org.apache.lucene.analysis.Tokenizer;
import org.apache.lucene.analysis.core.LetterTokenizer;
import org.apache.lucene.analysis.core.SimpleAnalyzer;
import org.apache.lucene.analysis.shingle.ShingleFilter;
import org.apache.lucene.analysis.tokenattributes.CharTermAttribute;
import org.apache.lucene.document.Document;
import org.apache.lucene.document.Field;
import org.apache.lucene.document.TextField;
import org.apache.lucene.index.DirectoryReader;
import org.apache.lucene.index.IndexReader;
import org.apache.lucene.index.IndexWriter;
import org.apache.lucene.index.IndexWriterConfig;
import org.apache.lucene.queryparser.classic.ParseException;
import org.apache.lucene.queryparser.classic.QueryParser;
import org.apache.lucene.search.IndexSearcher;
import org.apache.lucene.search.Query;
import org.apache.lucene.search.ScoreDoc;
import org.apache.lucene.search.spell.Dictionary;
import org.apache.lucene.search.spell.HighFrequencyDictionary;
import org.apache.lucene.search.spell.SpellChecker;
import org.apache.lucene.store.Directory;
import org.apache.lucene.store.RAMDirectory;
import org.jetbrains.annotations.NotNull;

public class TreatmentNameMatcher {
    private static final String DRUG_NAME_FIELD = "drugName";
    private static final String CANONICAL_DRUG_NAME_FIELD = "canonicalDrugName";
    private static final int NUM_HITS = 20;
    private static final int MAX_SHINGLES = 10;
    private static final float SPELLCHECK_ACCURACY = .8f;
    private static final float AMBIGUOUS_RESULTS_THRESHOLD = .7f;
    private final Directory index;
    private final SpellChecker spellChecker;

    public TreatmentNameMatcher(@NotNull final String mappingCsv) throws IOException {
        index = createIndex(mappingCsv);
        spellChecker = createIndexSpellchecker(index);
    }

    @NotNull
    public String match(@NotNull final String searchTerm) throws IOException, ParseException {
        final IndexReader reader = DirectoryReader.open(index);
        final IndexSearcher searcher = new IndexSearcher(reader);
        final Analyzer analyzer = spellcheckAnalyzer(spellChecker);
        final Query query = new QueryParser(DRUG_NAME_FIELD, analyzer).createPhraseQuery(DRUG_NAME_FIELD, searchTerm);
        final ScoreDoc[] hits = searcher.search(query, NUM_HITS).scoreDocs;
        if (hits.length > 0) {
            if (hits.length > 1) {
                final float topScoresSimilarity = hits[1].doc / hits[0].score;
                if (topScoresSimilarity > AMBIGUOUS_RESULTS_THRESHOLD) {
                    return "";
                }
            }
            return searcher.doc(hits[0].doc).get(CANONICAL_DRUG_NAME_FIELD);
        }
        return "";
    }

    @NotNull
    List<String> search(@NotNull final String searchTerm) throws IOException, ParseException {
        final String matchedTreatment = match(searchTerm);
        if (matchedTreatment.isEmpty()) {
            return matchMultiple(searchTerm);
        } else {
            return Lists.newArrayList(matchedTreatment);
        }
    }

    @NotNull
    private TokenStream getSpellcheckedShingleStream(@NotNull final String searchTerm) {
        StringReader reader = new StringReader(searchTerm);
        final Analyzer analyzer = createShingleAnalyzer(MAX_SHINGLES);
        return analyzer.tokenStream(DRUG_NAME_FIELD, reader);
    }

    @NotNull
    List<String> matchMultiple(@NotNull final String searchTerm) throws IOException, ParseException {
        final HashMap<String, String> tokenToTreatmentMap = Maps.newHashMap();
        final Set<String> matchedTokens = Sets.newHashSet();
        final TokenStream tokenStream = getSpellcheckedShingleStream(searchTerm);
        tokenStream.reset();
        while (tokenStream.incrementToken()) {
            final String searchToken = tokenStream.getAttribute(CharTermAttribute.class).toString();
            final String matchedTreatment = match(searchToken);
            if (!matchedTreatment.isEmpty()) {
                tokenToTreatmentMap.put(searchToken, matchedTreatment);
            }
        }
        tokenStream.end();
        tokenStream.close();
        tokenToTreatmentMap.keySet().stream().sorted(Comparator.comparing(String::length).reversed()).forEach(token -> {
            if (matchedTokens.stream().noneMatch(matchedToken -> matchedToken.contains(token))) {
                matchedTokens.add(token);
            }
        });
        return matchedTokens.stream().map(tokenToTreatmentMap::get).collect(Collectors.toList());
    }

    @NotNull
    private static Directory createIndex(@NotNull final String mappingCsv) throws IOException {
        final Directory treatmentIndex = new RAMDirectory();
        final CSVParser parser = CSVParser.parse(new File(mappingCsv), Charset.defaultCharset(), CSVFormat.DEFAULT.withHeader());
        final IndexWriter indexWriter = createIndexWriter(treatmentIndex);
        for (CSVRecord record : parser) {
            indexRecord(indexWriter, record);
        }
        indexWriter.close();
        return treatmentIndex;
    }

    @NotNull
    private static IndexWriter createIndexWriter(@NotNull final Directory directory) throws IOException {
        //        final Analyzer analyzer = createShingleAnalyzer(MAX_SHINGLES);
        final Analyzer analyzer = new SimpleAnalyzer();
        final IndexWriterConfig config = new IndexWriterConfig(analyzer);
        return new IndexWriter(directory, config);
    }

    @NotNull
    private static SpellChecker createIndexSpellchecker(@NotNull final Directory index) throws IOException {
        final Directory spellCheckerDirectory = new RAMDirectory();
        final IndexReader indexReader = DirectoryReader.open(index);
        final Analyzer analyzer = new SimpleAnalyzer();
        final IndexWriterConfig config = new IndexWriterConfig(analyzer);
        final Dictionary dictionary = new HighFrequencyDictionary(indexReader, DRUG_NAME_FIELD, 0.0f);
        final SpellChecker spellChecker = new SpellChecker(spellCheckerDirectory);
        spellChecker.indexDictionary(dictionary, config, false);
        spellChecker.setAccuracy(SPELLCHECK_ACCURACY);
        return spellChecker;
    }

    private static void indexRecord(@NotNull final IndexWriter writer, @NotNull final CSVRecord record) throws IOException {
        final String otherNamesString = record.get("other_names");
        final String canonicalName = record.get("drug");
        if (!otherNamesString.isEmpty()) {
            final CSVParser otherNamesParser = CSVParser.parse(otherNamesString, CSVFormat.DEFAULT);
            for (final CSVRecord otherNames : otherNamesParser) {
                for (final String name : otherNames) {
                    indexEntry(writer, name.trim(), canonicalName.trim());
                }
            }
        }
        indexEntry(writer, canonicalName.trim(), canonicalName.trim());
    }

    private static void indexEntry(@NotNull final IndexWriter writer, @NotNull final String name, @NotNull final String canonicalName)
            throws IOException {
        final Document document = new Document();
        document.add(new TextField(DRUG_NAME_FIELD, name, Field.Store.YES));
        document.add(new TextField(CANONICAL_DRUG_NAME_FIELD, canonicalName, Field.Store.YES));
        writer.addDocument(document);
    }

    @NotNull
    private static Analyzer createShingleAnalyzer(final int maxShingles) {
        return new Analyzer() {
            @Override
            protected TokenStreamComponents createComponents(final String field) {
                final StringReader reader = new StringReader(field);
                final Tokenizer source = new LetterTokenizer();
                source.setReader(reader);
                final TokenFilter filteredSource = new LowerCaseFilter(source);
                final ShingleFilter shingleFilter = new ShingleFilter(filteredSource, maxShingles);
                shingleFilter.setOutputUnigrams(true);
                return new TokenStreamComponents(source, shingleFilter);
            }
        };
    }

    @NotNull
    private static Analyzer spellcheckAnalyzer(@NotNull final SpellChecker spellChecker) {
        return new Analyzer() {
            @Override
            protected TokenStreamComponents createComponents(final String field) {
                final StringReader reader = new StringReader(field);
                final Tokenizer source = new LetterTokenizer();
                source.setReader(reader);
                final TokenFilter filteredSource = new LowerCaseFilter(source);
                final SpellCheckerTokenFilter spellCheckFilter = new SpellCheckerTokenFilter(filteredSource, spellChecker);
                return new TokenStreamComponents(source, spellCheckFilter);
            }
        };
    }
}
