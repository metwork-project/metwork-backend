REST_FRAMEWORK = {
	'DEFAULT_PERMISSION_CLASSES': (
	  'rest_framework.permissions.AllowAny',
	),
	'DEFAULT_AUTHENTICATION_CLASSES': (
		'rest_framework.authentication.TokenAuthentication',
		'rest_framework.authentication.SessionAuthentication',
	),
	'PAGE_SIZE': 100,
	'EXCEPTION_HANDLER': 'rest_framework_json_api.exceptions.exception_handler',
	'DEFAULT_PAGINATION_CLASS':
		'rest_framework_json_api.pagination.PageNumberPagination',
	'DEFAULT_PARSER_CLASSES': (
		'rest_framework_json_api.parsers.JSONParser',
		'rest_framework.parsers.FormParser',
		'rest_framework.parsers.MultiPartParser'
	),
	'DEFAULT_RENDERER_CLASSES': (
		'rest_framework_json_api.renderers.JSONRenderer',
		'rest_framework.renderers.BrowsableAPIRenderer',
	),
	'DEFAULT_METADATA_CLASS': 'rest_framework_json_api.metadata.JSONAPIMetadata',
}
