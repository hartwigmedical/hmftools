
plugins {
    id("com.hartwig.java-conventions")
    id("nu.studer.jooq")
    application
}

description = "HMF Tools - Patient Database"

dependencies {
    implementation(project(":hmf-common"))
    api(libs.jooq)
    jooqGenerator(libs.mysql.connector.java)
    implementation(libs.commons.dbcp2)
    implementation(libs.commons.csv)
    implementation(libs.lucene.core)
    implementation(libs.lucene.highlighter)
    implementation(libs.lucene.analyzers.common)
    implementation(libs.lucene.queryparser)
    implementation(libs.lucene.suggest)
    testImplementation(libs.junit)
    annotationProcessor(libs.immutables.value)
    compileOnly(libs.immutables.value)
}

// configure the jooq code generator
jooq {
    version.set(libs.jooq.get().versionConstraint.requiredVersion)

    configurations {
        create("main") {
            jooqConfiguration.apply {
                logging = org.jooq.meta.jaxb.Logging.WARN
                jdbc.apply {
                    driver = "com.mysql.cj.jdbc.Driver"
                    url = "jdbc:mysql://localhost:3306/hmfpatients_test?useSSL=false&serverTimezone=CET"
                    user = "build"
                    password = "build"
                }
                generator.apply {
                    name = "org.jooq.codegen.DefaultGenerator"
                    database.apply {
                        name = "org.jooq.meta.mysql.MySQLDatabase"
                        inputSchema = "hmfpatients_test"
                    }
                    target.apply {
                        packageName = "com.hartwig.hmftools.patientdb.database.hmfpatients"
                        directory = "${project.buildDir}/generated/sources/jooq"
                    }
                }
            }
        }
    }
}

application {
    mainClass.set("com.hartwig.hmftools.patientdb.clinical.LoadClinicalData")
}
