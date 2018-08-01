from base.models import User

if User.objects.count() == 0:

	users = [
		{
			'username': 'Guest',
			'email': 'metwork.dev@gmail.com',
			'password': 'AYL6jGBm6R',
		}
	]

	for ui in users:
		if User.objects.filter(username = ui['username']).count() == 0:
			u = User.objects.create(
				username = ui['username'],
				email = ui['email'])
			u.set_password(ui['password'])
			u.save()
