enableFeaturePreview("VERSION_CATALOGS")

rootProject.name = "hmftools"

include(":hmf-common")
include(":patient-db")
include(":gene-utils")
include(":stat-calcs")
include(":hmf-id-generator")
include(":fastq-stats")
include(":health-checker")
include(":paddle")
include(":compar")
include(":amber")
include(":cobalt")
include(":sage")
include(":gripss")
include(":purple")
include(":linx")
include(":sv-tools")
include(":virus-interpreter")
include(":lilac")
include(":neo")
include(":sigs")
include(":cuppa")
include(":isofox")
include(":pave")
include(":ckb-importer")
include(":iclusion-importer")
include(":vicc-importer")
include(":serve")
include(":protect")
include(":orange")
include(":patient-reporter")
include(":telo")

val javaVersion = settings.extra.get("java.version").toString()
val kotlinVersion = settings.extra.get("kotlin.version").toString()

pluginManagement {
    val kotlinVersion = settings.extra.get("kotlin.version").toString()
    plugins {
        kotlin("jvm") version kotlinVersion
        id("nu.studer.jooq") version "6.0.1"
    }
}

dependencyResolutionManagement {
    versionCatalogs {
        create("libs") {
            version("commons.cli", "1.3.1")
            version("immutables", "2.8.2")
            version("htsjdk", "2.23.0")
            version("intellij.annotations", "12.0")
            version("google.guava", "31.0.1-jre")
            version("google.gson", "2.8.9")
            version("apache.commons.lang3", "3.6")
            version("apache.commons.math3", "3.6")
            version("apache.commons.csv", "1.5")
            version("apache.commons.dbcp2", "2.1.1")
            version("apache.log4j", "2.17.1")
            version("apache.lucene", "7.1.0")
            version("jooq", "3.15.5")
            version("mysqlconnector", "8.0.16")
            version("rxjava2", "2.1.2")
            version("retrofit", "2.4.0")
            version("moshi", "1.6.0")
            version("okhttp", "3.9.0")
            version("bouncycastle.jdk15", "1.53")
            version("kotlin", kotlinVersion)
            version("kotlin.coroutines", "1.5.2")
            version("selenium", "3.14.0")
            version("itext", "7.1.5")
            version("tablesaw", "0.38.5")
            version("jcommander", "1.81")

            version("kotest", "4.6.4")
            version("junit", "4.13.1")
            version("mockito", "2.23.4")

            alias("htsjdk").to("com.github.samtools", "htsjdk").versionRef("htsjdk")
            alias("annotations").to("com.intellij", "annotations").versionRef("intellij.annotations")
            alias("google.guava").to("com.google.guava", "guava").versionRef("google.guava")
            alias("google.gson").to("com.google.code.gson", "gson").versionRef("google.gson")
            alias("commons.cli").to("commons-cli", "commons-cli").versionRef("commons.cli")
            alias("commons.lang3").to("org.apache.commons", "commons-lang3").versionRef("apache.commons.lang3")
            alias("commons.math3").to("org.apache.commons", "commons-math3").versionRef("apache.commons.math3")
            alias("commons.csv").to("org.apache.commons", "commons-csv").versionRef("apache.commons.csv")
            alias("commons.dbcp2").to("org.apache.commons", "commons-dbcp2").versionRef("apache.commons.dbcp2")
            alias("log4j.core").to("org.apache.logging.log4j", "log4j-core").versionRef("apache.log4j")
            alias("log4j.slf4j.impl").to("org.apache.logging.log4j", "log4j-slf4j-impl").versionRef("apache.log4j")
            alias("jooq").to("org.jooq", "jooq").versionRef("jooq")
            alias("mysql.connector.java").to("mysql", "mysql-connector-java").versionRef("mysqlconnector")
            alias("immutables.value").to("org.immutables", "value").versionRef("immutables")
            alias("immutables.gson").to("org.immutables", "gson").versionRef("immutables")
            alias("rxjava").to("io.reactivex.rxjava2", "rxjava").versionRef("rxjava2")
            alias("okhttp").to("com.squareup.okhttp3", "okhttp").versionRef("okhttp")
            alias("retrofit").to("com.squareup.retrofit2", "retrofit").versionRef("retrofit")
            alias("adapter.rxjava2").to("com.squareup.retrofit2", "adapter-rxjava2").versionRef("retrofit")
            alias("converter.gson").to("com.squareup.retrofit2", "converter-gson").versionRef("retrofit")
            alias("converter.moshi").to("com.squareup.retrofit2", "converter-moshi").versionRef("retrofit")
            alias("moshi.core").to("com.squareup.moshi", "moshi").versionRef("moshi")
            alias("moshi.kotlin").to("com.squareup.moshi", "moshi-kotlin").versionRef("moshi")
            alias("lucene.core").to("org.apache.lucene", "lucene-core").versionRef("apache.lucene")
            alias("lucene.analyzers-common").to("org.apache.lucene", "lucene-analyzers-common").versionRef("apache.lucene")
            alias("lucene.queryparser").to("org.apache.lucene", "lucene-queryparser").versionRef("apache.lucene")
            alias("lucene.suggest").to("org.apache.lucene", "lucene-suggest").versionRef("apache.lucene")
            alias("lucene.highlighter").to("org.apache.lucene", "lucene-highlighter").versionRef("apache.lucene")
            alias("bcprov.jdk15on").to("org.bouncycastle", "bcprov-jdk15on").versionRef("bouncycastle.jdk15")
            alias("kotlin.reflect").to("org.jetbrains.kotlin", "kotlin-reflect").versionRef("kotlin")
            alias("kotlinx.coroutines.core").to("org.jetbrains.kotlinx", "kotlinx-coroutines-core").versionRef("kotlin.coroutines")
            alias("selenium.java").to("org.seleniumhq.selenium", "selenium-java").versionRef("selenium")
            alias("kernel").to("com.itextpdf", "kernel").versionRef("itext")
            alias("layout").to("com.itextpdf", "layout").versionRef("itext")
            alias("io").to("com.itextpdf", "io").versionRef("itext")
            alias("tablesaw.core").to("tech.tablesaw", "tablesaw-core").versionRef("tablesaw")
            alias("jcommander").to("com.beust", "jcommander").versionRef("jcommander")
            alias("mockito.core").to("org.mockito", "mockito-core").versionRef("mockito")
            alias("junit").to("junit", "junit").versionRef("junit")
            alias("kotest.runner.junit5.jvm").to("io.kotest", "kotest-runner-junit5-jvm").versionRef("kotest")
            alias("kotest.assertions.core.jvm").to("io.kotest", "kotest-assertions-core-jvm").versionRef("kotest")
        }
    }
}

// add a Maven style build summary report
plugins {
    id("org.kordamp.gradle.insight") version "0.47.0"
}