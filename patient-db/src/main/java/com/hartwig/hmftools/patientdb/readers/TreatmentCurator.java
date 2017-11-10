package com.hartwig.hmftools.patientdb.readers;

import java.io.File;
import java.io.IOException;
import java.io.StringReader;
import java.nio.charset.Charset;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.patientdb.data.CuratedTreatment;
import com.hartwig.hmftools.patientdb.data.ImmutableCuratedTreatment;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.lucene.analysis.Analyzer;
import org.apache.lucene.analysis.LowerCaseFilter;
import org.apache.lucene.analysis.TokenFilter;
import org.apache.lucene.analysis.TokenStream;
import org.apache.lucene.analysis.Tokenizer;
import org.apache.lucene.analysis.core.SimpleAnalyzer;
import org.apache.lucene.analysis.shingle.ShingleFilter;
import org.apache.lucene.analysis.tokenattributes.CharTermAttribute;
import org.apache.lucene.analysis.util.CharTokenizer;
import org.apache.lucene.document.Document;
import org.apache.lucene.document.Field;
import org.apache.lucene.document.StringField;
import org.apache.lucene.document.TextField;
import org.apache.lucene.index.DirectoryReader;
import org.apache.lucene.index.IndexReader;
import org.apache.lucene.index.IndexWriter;
import org.apache.lucene.index.IndexWriterConfig;
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

public class TreatmentCurator {
    private static final Logger LOGGER = LogManager.getLogger(TreatmentCurator.class);
    private static final String DRUG_NAME_FIELD = "drugName";
    private static final String CANONICAL_DRUG_NAME_FIELD = "canonicalDrugName";
    private static final String DRUG_TYPE_FIELD = "drugType";
    private static final int NUM_HITS = 20;
    private static final int MAX_SHINGLES = 10;
    private static final float SPELLCHECK_ACCURACY = .75f;
    private static final float AMBIGUOUS_RESULTS_THRESHOLD = .7f;
    private final SpellChecker spellChecker;
    private final IndexSearcher indexSearcher;

    public TreatmentCurator(@NotNull final String mappingCsv) throws IOException {
        final Directory index = createIndex(mappingCsv);
        final IndexReader reader = DirectoryReader.open(index);
        spellChecker = createIndexSpellchecker(index);
        indexSearcher = new IndexSearcher(reader);
    }

    @NotNull
    public List<CuratedTreatment> search(@NotNull final String searchTerm) throws IOException {
        final Optional<CuratedTreatment> matchedTreatment = matchSingle(searchTerm);
        if (!matchedTreatment.isPresent()) {
            LOGGER.warn("Failed to match search term: {}. attempting to match multiple treatments", searchTerm);
            final List<CuratedTreatment> matchedTreatments = matchMultiple(searchTerm);
            if (!matchedTreatments.isEmpty()) {
                LOGGER.info("Matched multiple treatments {} to {}", searchTerm,
                        matchedTreatments.stream().map(CuratedTreatment::name).collect(Collectors.toList()));
            }
            return matchMultiple(searchTerm);
        } else {
            return Lists.newArrayList(matchedTreatment.get());
        }
    }

    @NotNull
    Optional<CuratedTreatment> matchSingle(@NotNull final String searchTerm) throws IOException {
        final Analyzer analyzer = spellcheckAnalyzer(spellChecker);
        final Query query = new QueryParser(DRUG_NAME_FIELD, analyzer).createPhraseQuery(DRUG_NAME_FIELD, searchTerm);
        final ScoreDoc[] hits = indexSearcher.search(query, NUM_HITS).scoreDocs;
        if (hits.length > 0) {
            return filterSearchResults(hits);
        }
        return Optional.empty();
    }

    @NotNull
    private Optional<CuratedTreatment> filterSearchResults(@NotNull final ScoreDoc[] hits) throws IOException {
        final ScoreDoc topScoreHit = hits[0];
        final Document topSearchResult = indexSearcher.doc(hits[0].doc);
        final String topResultCanonicalDrug = topSearchResult.get(CANONICAL_DRUG_NAME_FIELD);
        for (final ScoreDoc hit : hits) {
            final String hitCanonicalDrug = indexSearcher.doc(hit.doc).get(CANONICAL_DRUG_NAME_FIELD);
            if (!hitCanonicalDrug.equals(topResultCanonicalDrug)) {
                final float scoreSimilarity = hit.score / topScoreHit.score;
                if (scoreSimilarity > AMBIGUOUS_RESULTS_THRESHOLD) {
                    LOGGER.info("Search results ambiguous. topHit: {}, otherHit: {}", topResultCanonicalDrug,
                            indexSearcher.doc(hit.doc).get(CANONICAL_DRUG_NAME_FIELD));
                    return Optional.empty();
                }
            }
        }
        return Optional.of(ImmutableCuratedTreatment.of(topResultCanonicalDrug, topSearchResult.get(DRUG_TYPE_FIELD)));
    }

    @NotNull
    List<CuratedTreatment> matchMultiple(@NotNull final String searchTerm) throws IOException {
        final HashMap<String, CuratedTreatment> tokenToTreatmentMap = Maps.newHashMap();
        final Set<String> matchedTokens = Sets.newHashSet();
        final TokenStream tokenStream = getSpellcheckedShingleStream(searchTerm);
        tokenStream.reset();
        while (tokenStream.incrementToken()) {
            final String searchToken = tokenStream.getAttribute(CharTermAttribute.class).toString();
            final Optional<CuratedTreatment> matchedTreatment = matchSingle(searchToken);
            matchedTreatment.ifPresent(curatedTreatment -> tokenToTreatmentMap.put(searchToken, curatedTreatment));
        }
        tokenStream.end();
        tokenStream.close();
        tokenToTreatmentMap.keySet().stream().sorted(Comparator.comparing(String::length).reversed()).forEach(token -> {
            if (matchedTokens.stream().noneMatch(matchedToken -> matchedToken.contains(token))) {
                matchedTokens.add(token);
            }
        });
        return matchedTokens.stream().map(tokenToTreatmentMap::get).distinct().collect(Collectors.toList());
    }

    @NotNull
    private TokenStream getSpellcheckedShingleStream(@NotNull final String searchTerm) {
        StringReader reader = new StringReader(searchTerm);
        final Analyzer analyzer = createShingleAnalyzer(MAX_SHINGLES);
        return analyzer.tokenStream(DRUG_NAME_FIELD, reader);
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
        final Analyzer analyzer = createIndexAnalyzer();
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
        final String drugType = record.get("type");
        final String canonicalName = record.get("drug");
        if (!otherNamesString.isEmpty()) {
            final CSVParser otherNamesParser = CSVParser.parse(otherNamesString, CSVFormat.DEFAULT);
            for (final CSVRecord otherNames : otherNamesParser) {
                for (final String name : otherNames) {
                    indexEntry(writer, name.trim(), drugType, canonicalName.trim());
                }
            }
        }
        indexEntry(writer, canonicalName.trim(), drugType, canonicalName.trim());
    }

    private static void indexEntry(@NotNull final IndexWriter writer, @NotNull final String name, @NotNull final String type,
            @NotNull final String canonicalName) throws IOException {
        final Document document = new Document();
        document.add(new TextField(DRUG_NAME_FIELD, name, Field.Store.YES));
        document.add(new StringField(DRUG_TYPE_FIELD, type, Field.Store.YES));
        document.add(new TextField(CANONICAL_DRUG_NAME_FIELD, canonicalName, Field.Store.YES));
        writer.addDocument(document);
    }

    @NotNull
    private static Analyzer createShingleAnalyzer(final int maxShingles) {
        return new Analyzer() {
            @Override
            protected TokenStreamComponents createComponents(final String field) {
                final StringReader reader = new StringReader(field);
                final Tokenizer source = CharTokenizer.fromTokenCharPredicate(Character::isLetterOrDigit);
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
                final Tokenizer source = CharTokenizer.fromTokenCharPredicate(Character::isLetterOrDigit);
                source.setReader(reader);
                final TokenFilter filteredSource = new LowerCaseFilter(source);
                final SpellCheckerTokenFilter spellCheckFilter = new SpellCheckerTokenFilter(filteredSource, spellChecker);
                return new TokenStreamComponents(source, spellCheckFilter);
            }
        };
    }

    @NotNull
    private static Analyzer createIndexAnalyzer() {
        return new Analyzer() {
            @Override
            protected TokenStreamComponents createComponents(final String field) {
                final StringReader reader = new StringReader(field);
                final Tokenizer source = CharTokenizer.fromTokenCharPredicate(Character::isLetterOrDigit);
                source.setReader(reader);
                final TokenFilter filteredSource = new LowerCaseFilter(source);
                return new TokenStreamComponents(source, filteredSource);
            }
        };
    }
}
