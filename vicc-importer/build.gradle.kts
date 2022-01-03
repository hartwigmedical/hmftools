
plugins {
    id("com.hartwig.java-conventions")
    id("nu.studer.jooq")
    application
}

description = "HMF Tools - VICC Importer"

dependencies {
    implementation(project(":hmf-common"))
    implementation(libs.jooq)
    jooqGenerator(libs.mysql.connector.java)
    annotationProcessor(libs.immutables.value)
    compileOnly(libs.immutables.value)
    testImplementation(libs.junit)
}

// configure the jooq generator
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
                        inputSchema = "vicc_test"
                    }
                    target.apply {
                        packageName = "com.hartwig.hmftools.vicc.database"
                        directory = "${project.buildDir}/generated/sources/jooq"
                    }
                }
            }
        }
    }
}

application {
    mainClass.set("com.hartwig.hmftools.vicc.ViccJsonSQLImporter")
}
