
plugins {
    id("com.hartwig.java-conventions")
    id("nu.studer.jooq")
}

description = "HMF Tools - CKB Importer"

dependencies {
    implementation(project(":hmf-common"))
    implementation(libs.jooq)
    jooqGenerator(libs.mysql.connector.java)
    annotationProcessor(libs.immutables.value)
    compileOnly(libs.immutables.value)
    testImplementation(libs.junit)
}

shadowJar {
    mainClass.set("com.hartwig.hmftools.ckb.CkbImporterApplication")
}

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
                        inputSchema = "ckb_test"
                    }
                    target.apply {
                        packageName = "com.hartwig.hmftools.ckb.database"
                        directory = "${project.buildDir}/generated/sources/jooq"
                    }
                }
            }
        }
    }
}

