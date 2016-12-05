import soundDownload as SD

key = 'CMaNlgqF7Ba6N39NZhnQcQE0Dn7DlN1mJWo9GRIO'

SD.downloadSoundsFreesound(queryText='flute-G6', tag='single-note', duration=(0,10), API_Key=key, outputDir='tmp', topNResults=1, featureExt='.json')

SD.downloadSoundsFreesound(queryText='cello-B2', tag='single-note', duration=(0,10), API_Key=key, outputDir='tmp', topNResults=1, featureExt='.json')

SD.downloadSoundsFreesound(queryText='Guitar_E', tag='single-note', duration=(0,10), API_Key=key, outputDir='tmp', topNResults=1, featureExt='.json')
